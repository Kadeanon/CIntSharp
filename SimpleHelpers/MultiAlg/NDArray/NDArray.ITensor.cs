using NetFabric.Numerics.Tensors.Operators;
using SimpleHelpers.MultiAlg.Helpers;
using SimpleHelpers.MultiAlg.SimpleTranspose;
using SimpleHelpers.Utilities;
using System.Buffers;
using System.Collections;
using System.Numerics.Tensors;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray : ITensor<NDArray, double>
    {
        #region Many indexer
        object? ITensor.
            this[params ReadOnlySpan<nint> indexes]
        {
            get =>
                indexes.Length == Rank ?
                this[indexes] : Slice(indexes);

            set
            {
                if (indexes.Length == Rank)
                {
                    if(value is IConvertible conv)
                        this[indexes] = 
                            conv.ToDouble(Thread.CurrentThread.CurrentCulture);
                    else
                        throw new ArgumentException(
                            "The value must be convertible to double.",
                            nameof(value));
                }
                else
                {
                    if (value is NDArray array)
                    {
                        TransMethods.CopyTo(array, Slice(indexes));
                    }
                    else
                    {
                        throw new ArgumentException(
                            "The value must be an NDArray or convertible to double.",
                            nameof(value));
                    }
                }
            }
        }

        object? ITensor.
            this[params ReadOnlySpan<LengthIndex> indexes]
        {
            get =>
                indexes.Length == Rank ?
                this[indexes] : Slice(indexes);

            set
            {
                if (indexes.Length == Rank)
                {
                    if (value is IConvertible conv)
                        this[indexes] =
                            conv.ToDouble(Thread.CurrentThread.CurrentCulture);
                    else
                        throw new ArgumentException(
                            "The value must be convertible to double.",
                            nameof(value));
                }
                else
                {
                    if (value is NDArray array)
                    {
                        TransMethods.CopyTo(array, Slice(indexes));
                    }
                    else
                    {
                        throw new ArgumentException(
                            "The value must be an NDArray or convertible to double.",
                            nameof(value));
                    }
                }
            }
        }

        ref readonly double IReadOnlyTensor<NDArray, double>.
            this[params ReadOnlySpan<nint> indexes] =>
            ref this[indexes];

        ref readonly double IReadOnlyTensor<NDArray, double>.
            this[params ReadOnlySpan<LengthIndex> indexes] =>
            ref this[indexes];

        object? IReadOnlyTensor.
            this[params ReadOnlySpan<nint> indexes] =>
                indexes.Length == Rank ?
                this[indexes] : Slice(indexes);

        object? IReadOnlyTensor.
            this[params ReadOnlySpan<LengthIndex> indexes] =>
                indexes.Length == Rank ?
                this[indexes] : Slice(indexes);
        #endregion Many indexer

        #region Properties
        public static NDArray Empty { get; } = new([], []);

        public bool IsEmpty => Size == 0;

        public bool IsReadOnly => false;

        nint IReadOnlyTensor.FlattenedLength => Size;

        public bool HasAnyDenseDimensions =>
            Rank != Continuous.NumLayer;

        public bool IsDense => !State.HasFlag(ArrayState.Broken);

        public bool IsPinned => false;
        #endregion Properties

        static NDArray ITensor<NDArray, double>.Create(
            ReadOnlySpan<nint> lengths, bool pinned)
        {
            if (pinned)
                throw new NotSupportedException(
                    "Pinned NDArray creation is not supported.");
            return Create(lengths);
        }

        static NDArray ITensor<NDArray, double>.Create(
            ReadOnlySpan<nint> lengths, ReadOnlySpan<nint> strides,
            bool pinned)
        {
            if (pinned)
                throw new NotSupportedException(
                    "Pinned NDArray creation is not supported.");
            return Create(lengths, strides);
        }

        static NDArray ITensor<NDArray, double>.CreateUninitialized
            (ReadOnlySpan<nint> lengths, bool pinned)
        {
            if (pinned)
                throw new NotSupportedException(
                    "Pinned NDArray creation is not supported.");
            return CreateUninitialized(lengths);
        }

        static NDArray ITensor<NDArray, double>.CreateUninitialized
            (ReadOnlySpan<nint> lengths, ReadOnlySpan<nint> strides,
            bool pinned)
        {
            if (pinned)
                throw new NotSupportedException(
                    "Pinned NDArray creation is not supported.");
            return CreateUninitialized(lengths, strides);
        }

        public ReadOnlyTensorSpan<double> AsReadOnlyTensorSpan()
        {
            return new ReadOnlyTensorSpan<double>(
                Data, (int)Offset, Lengths, Strides);
        }

        public ReadOnlyTensorSpan<double> AsReadOnlyTensorSpan(
            params ReadOnlySpan<nint> startIndexes)
        {
            throw new NotImplementedException();
        }

        public ReadOnlyTensorSpan<double> AsReadOnlyTensorSpan(
            params ReadOnlySpan<LengthIndex> startIndexes)
        {
            throw new NotImplementedException();
        }

        public ReadOnlyTensorSpan<double> AsReadOnlyTensorSpan(
            params ReadOnlySpan<NRange> ranges)
        {
            throw new NotImplementedException();
        }

        public TensorSpan<double> AsTensorSpan()
        {
            return new (
                Data, (int)Offset, Lengths, Strides);
        }

        public TensorSpan<double> AsTensorSpan(params ReadOnlySpan<nint> startIndexes)
        {
            throw new NotImplementedException();
        }

        public TensorSpan<double> AsTensorSpan(params ReadOnlySpan<LengthIndex> startIndexes)
        {
            throw new NotImplementedException();
        }

        public TensorSpan<double> AsTensorSpan(params ReadOnlySpan<NRange> ranges)
        {
            return Slice(ranges).AsTensorSpan();
        }

        public void Clear()
            => ApplySelf<IdentityOperator<double>, double>(this, 0.0);

        public void CopyTo(in TensorSpan<double> destination)
        {
            throw new NotImplementedException();
        }

        public void Fill(double value)
            => ApplySelf<IdentityOperator<double>, double>(this, value);

        public void Fill(object value)
        {
            if (value is IConvertible conv)
                Fill(conv.ToDouble(Thread.CurrentThread.CurrentCulture));
            else
                throw new ArgumentException(
                    "The value must be convertible to double.",
                    nameof(value));
        }

        public void FlattenTo(Span<double> destination)
        {
            AsReadOnlyTensorSpan().FlattenTo(destination);
        }

        public IEnumerator<double> GetEnumerator()
        {
            throw new NotImplementedException();
        }

        public ref double GetPinnableReference()
        {
            return ref GetHeadRef();
        }

        public MemoryHandle GetPinnedHandle()
        {
            throw new NotImplementedException();
        }

        public NDArray Slice(params ReadOnlySpan<nint> startIndexes)
        {
            ArgumentOutOfRangeException.ThrowIfGreaterThan(
                startIndexes.Length, Rank,
                $"The number of start indexes ({startIndexes.Length}) " +
                $"must be less or equal to the rank of the array ({Rank}).");
            nint offset = Offset;
            nint[] metadata = new nint[Rank * 2];
            int i = 0;
            for (; i < startIndexes.Length; i++)
            {
                var dimIndex = startIndexes[i];
                var dimTotalLength = Lengths[i];
                ThrowUtils.ThrowIfNotInRange_CheckNegative(ref dimIndex,
                     dimTotalLength, nameof(startIndexes));
                offset += dimIndex * Strides[i];
                metadata[i] = dimTotalLength - dimIndex;
                metadata[Rank + i] = Strides[i];
            }
            for (; i < Rank; i++)
            {
                metadata[i] = Lengths[i];
                metadata[Rank + i] = Strides[i];
            }
            return new NDArray(Data, offset, metadata, Rank);
        }

        public NDArray Slice(params ReadOnlySpan<LengthIndex> startIndexes)
        {
            ArgumentOutOfRangeException.ThrowIfGreaterThan(
                startIndexes.Length, Rank,
                $"The number of start indexes ({startIndexes.Length}) " +
                $"must be less or equal to the rank of the array ({Rank}).");
            nint offset = Offset;
            nint[] metadata = new nint[Rank * 2];
            int i = 0;
            for (; i < startIndexes.Length; i++)
            {
                var dimIndex = startIndexes[i].GetOffset(Lengths[i]);
                var dimTotalLength = Lengths[i];
                ThrowUtils.ThrowIfNotInRange_LeftClosedRightOpen(dimIndex,
                    0, dimTotalLength, nameof(startIndexes));
                offset += dimIndex * Strides[i];
                metadata[i] = dimTotalLength - dimIndex;
                metadata[Rank + i] = Strides[i];
            }
            for (; i < Rank; i++)
            {
                metadata[i] = Lengths[i];
                metadata[Rank + i] = Strides[i];
            }
            return new NDArray(Data, offset, metadata, Rank);
        }

        public NDArray Slice(params ReadOnlySpan<NRange> ranges)
        {
            ArgumentOutOfRangeException.ThrowIfGreaterThan(
                ranges.Length, Rank,
                $"The number of ranges ({ranges.Length}) " +
                $"must be less or equal to the rank of the array ({Rank}).");
            nint offset = Offset;
            nint[] metadata = new nint[Rank * 2];
            int i = 0;
            for (; i < ranges.Length; i++)
            {
                var range = ranges[i];
                var dimTotalLength = Lengths[i];
                var dimStride = Strides[i];
                var (dimStart, dimLength) = ThrowUtils.ThrowIfNotInRange_CheckNRange(
                        range, dimTotalLength, nameof(ranges));
                offset += dimStart * dimStride;
                metadata[i] = dimLength;
                metadata[Rank + i] = dimStride;
            }
            for (; i < Rank; i++)
            {
                metadata[i] = Lengths[i];
                metadata[Rank + i] = Strides[i];
            }
            return new NDArray(Data, offset, metadata, Rank);
        }

        public NDArray ToDenseTensor()
        {
            if (IsDense)
                return this;
            return Reshape(Lengths, forceCopy:true);
        }

        public bool TryCopyTo(in TensorSpan<double> destination)
        {
            throw new NotImplementedException();
        }

        public bool TryFlattenTo(Span<double> destination)
        {
            throw new NotImplementedException();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            throw new NotImplementedException();
        }

        ref readonly double IReadOnlyTensor<NDArray, double>.GetPinnableReference()
        {
            throw new NotImplementedException();
        }
    }
}
