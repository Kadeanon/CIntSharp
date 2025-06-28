using SimpleHelpers.MultiAlg.Helpers;
using SimpleHelpers.Utilities.Pools;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Text;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray 
    {
        #region Menbers
        internal readonly double[] Data;

        internal nint Offset { get; }
        #endregion

        #region Properties

        internal nint[] Metadata { get; }

        public int Rank { get; }

        public ContinuousInfo Continuous { get; }

        /// <summary>
        /// Indecate the state of the array.
        /// </summary>
        /// <remarks>
        /// It is a combination of <see cref="ArrayState"/> flags. So use <see cref="Enum.HasFlag(Enum)"/> to check the state.
        /// </remarks>
        public ArrayState State
        {
            get
            {
                ArrayState state = ArrayState.Default;
                ContinuousInfo info = Continuous;
                if (info.Layers[0].Stride == 1)
                    state |= ArrayState.FortranStyle;

                if (info.Layers[^1].Stride != 1)
                {
                    state |= ArrayState.Broken;
                    state |= ArrayState.Segmented;
                }
                else
                {
                    state |= ArrayState.CStyle;
                    if (info.NumLayer != 1)
                    {
                        state |= ArrayState.Segmented;
                    }
                }
                return state;
            }
        }

        public ref double GetHeadRef()
        {
            if(Rank == 0)
                return ref Unsafe.NullRef<double>();
            return ref Data[Offset];
        }

        public ReadOnlySpan<nint> Lengths => Metadata.AsSpan(0, Rank);

        public ReadOnlySpan<nint> Strides => Metadata.AsSpan(Rank, Rank);

        public nint Size => TensorPrimitives.Product(Lengths);

        public bool IsScalar => Size == 1;
        #endregion

        #region Constructors
        public NDArray(double[] data, ReadOnlySpan<nint> shape)
        {
            Data = data;
            Offset = 0;
            Rank = shape.Length;
            var metadata = new nint[Rank * 2];
            Metadata = metadata;
            if (shape.Length == 0)
            {
                Continuous = new(metadata);
                return;
            }
            nint stride = 1;
            for (int i = shape.Length - 1; i >= 0; i--)
            {
                nint length = shape[i];
                metadata[i] = length;
                metadata[Rank + i] = stride;
                stride *= length;
            }
            ArgumentOutOfRangeException.ThrowIfGreaterThan(
                stride, data.Length, nameof(shape));
            Continuous = new(metadata);
        }

        public NDArray(double[] data,
            nint offset, ReadOnlySpan<nint> shape)
        {
            Data = data;
            Offset = offset;
            Rank = shape.Length;
            var metadata = new nint[Rank * 2];
            Metadata = metadata;
            if (shape.Length == 0)
            {
                Continuous = new(metadata);
                return;
            }
            nint stride = 1;
            for (int i = shape.Length - 1; i >= 0; i--)
            {
                nint length = shape[i];
                metadata[i] = length;
                metadata[Rank + i] = stride;
                stride *= length;
            }
            ArgumentOutOfRangeException.ThrowIfGreaterThan(
                offset + stride, data.Length, nameof(shape));
            Continuous = new(metadata);
        }

        public NDArray(double[] data,
            ReadOnlySpan<nint> shape, ReadOnlySpan<nint> strides)
        {
            Data = data;
            Offset = 0;
            Rank = shape.Length;
            var metadata = new nint[Rank * 2];
            Metadata = metadata;
            if (shape.Length == 0)
            {
                Continuous = new(metadata);
                return;
            }
            nint totalLength = 1;
            for (int i = shape.Length - 1; i >= 0; i--)
            {
                nint length = shape[i];
                nint stride = strides[i];
                metadata[i] = length;
                metadata[Rank + i] = stride;
                totalLength += (length - 1) * stride;
            }
            ArgumentOutOfRangeException.ThrowIfGreaterThan(
                totalLength, data.Length, nameof(shape));
            Continuous = new(metadata);
        }

        public NDArray(double[] data,
            nint offset, ReadOnlySpan<nint> shape, ReadOnlySpan<nint> strides)
        {
            Data = data;
            Offset = offset;
            Rank = shape.Length;
            var metadata = new nint[Rank * 2];
            Metadata = metadata;
            if (shape.Length == 0)
            {
                Continuous = new(metadata);
                return;
            }
            nint totalLength = 1;
            for (int i = shape.Length - 1; i >= 0; i--)
            {
                nint length = shape[i];
                nint stride = strides[i];
                metadata[i] = length;
                metadata[Rank + i] = stride;
                totalLength += (length - 1) * stride;
            }
            ArgumentOutOfRangeException.ThrowIfGreaterThan(
                offset + totalLength, data.Length, nameof(shape));
            Continuous = new(metadata);
        }

        public NDArray(double[] data, nint offset, nint[] metadata, int rank)
        {
            Data = data;
            Offset = offset;
            Metadata = metadata;
            Rank = rank;
            Continuous = new(metadata);
        }
        #endregion Constructors

        #region Alloc
        public static NDArray CreateUninitialized(ReadOnlySpan<nint> shape)
        {
            nint length = TensorPrimitives.Product(shape);
            if (length > Array.MaxLength)
                throw new ArgumentOutOfRangeException(nameof(shape),
                    "The array is too large!");
            double[] data = GC.AllocateUninitializedArray<double>((int)length);
            return new NDArray(data, shape);
        }

        public static NDArray Create(params ReadOnlySpan<nint> shape)
        {
            nint length = TensorPrimitives.Product(shape);
            double[] data = new double[length];
            return new NDArray(data, shape);
        }

        public static NDArray Create(ReadOnlySpan<nint> shape, double val)
        {
            nint length = TensorPrimitives.Product(shape);
            double[] data = new double[length];
            data.AsSpan().Fill(val);
            return new NDArray(data, shape);
        }

        public static NDArray CreateUninitialized(ReadOnlySpan<nint> shape,
            ReadOnlySpan<nint> strides)
        {
            ArgumentOutOfRangeException.ThrowIfNotEqual
                (shape.Length, strides.Length);
            int rank = shape.Length;
            var metadata = new nint[rank * 2];
            nint totalLength = 1;
            for (int i = 0; i < shape.Length; i++)
            {
                nint length = shape[i];
                nint stride = strides[i];
                metadata[i] = length;
                metadata[rank + i] = stride;
                totalLength += (length - 1) * stride;
            }
            if (totalLength > Array.MaxLength)
                throw new ArgumentOutOfRangeException(nameof(shape),
                    "The array is too large!");
            double[] data = GC.
                AllocateUninitializedArray<double>((int)totalLength);
            return new NDArray(data, 0, metadata, rank);
        }

        public static NDArray Create(ReadOnlySpan<nint> shape,
            ReadOnlySpan<nint> strides)
        {
            ArgumentOutOfRangeException.ThrowIfNotEqual
                (shape.Length, strides.Length);
            int rank = shape.Length;
            var metadata = new nint[rank * 2];
            nint totalLength = 1;
            for (int i = 0; i < shape.Length; i++)
            {
                nint length = shape[i];
                nint stride = strides[i];
                metadata[i] = length;
                metadata[rank + i] = stride;
                totalLength += (length - 1) * stride;
            }
            if (totalLength > Array.MaxLength)
                throw new ArgumentOutOfRangeException(nameof(shape),
                    "The array is too large!");
            double[] data = new double[totalLength];
            return new NDArray(data, 0, metadata, rank);
        }

        public static NDArray Create(ReadOnlySpan<nint> shape,
            ReadOnlySpan<nint> strides, double val)
        {
            ArgumentOutOfRangeException.ThrowIfNotEqual
                (shape.Length, strides.Length);
            int rank = shape.Length;
            var metadata = new nint[rank * 2];
            nint totalLength = 1;
            for (int i = 0; i < shape.Length; i++)
            {
                nint length = shape[i];
                nint stride = strides[i];
                metadata[i] = length;
                metadata[rank + i] = stride;
                totalLength += (length - 1) * stride;
            }
            if (totalLength > Array.MaxLength)
                throw new ArgumentOutOfRangeException(nameof(shape),
                    "The array is too large!");
            double[] data = new double[totalLength];
            data.AsSpan().Fill(val);
            return new NDArray(data, 0, metadata, rank);
        }

        public NDArray UninitializedLike()
            => CreateUninitialized(Lengths);

        public static NDArray UninitializedLike(NDArray array)
            => CreateUninitialized(array.Lengths);

        public NDArray ZeroLike()
            => Create(Lengths);

        public static NDArray ZeroLike(NDArray array)
            => Create(array.Lengths);

        public NDArray FillLike(double val)
            => Create(Lengths, val);

        public static NDArray FillLike(NDArray array, double val)
            => Create(array.Lengths, val);

        public NDArray Clone()
        {
            var copied = UninitializedLike();
            SimpleTranspose.TransMethods.CopyTo(this, copied);
            return copied;
        }
        #endregion Alloc


        #region String

        public string MetaDataString()
        {
            using (StringBuilderPool.Borrow(out var sb))
            {
                MetaDataString(sb);
                return sb.ToString();
            }
        }

        public void MetaDataString(StringBuilder sb)
        {
            sb
            .AppendLine($"Shape: ({string.Join(", ", Lengths.ToArray())})")
            .AppendLine($"Data Type: {typeof(double).Name}")
            .AppendLine($"Memory Usage: {Size * sizeof(double)} Byte")
            ;
        }
        #endregion

    }
}
