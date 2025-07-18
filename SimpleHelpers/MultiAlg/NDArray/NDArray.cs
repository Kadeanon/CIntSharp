using SimpleHelpers.LinearAlg;
using SimpleHelpers.MultiAlg.Helpers;
using SimpleHelpers.Utilities;
using SimpleHelpers.Utilities.Pools;
using System.Buffers;
using System.Diagnostics.CodeAnalysis;
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
            long totalSize = length * 8;
            var sizeB = totalSize % 1024;
            totalSize /= 1024;
            var sizeKB = totalSize % 1024;
            totalSize /= 1024;
            var sizeMB = totalSize % 1024;
            totalSize /= 1024;
            var sizeGB = totalSize;
            if (length * 8 > Array.MaxLength)
            {
                throw new ArgumentOutOfRangeException(nameof(shape),
                    "The array is too large! It need to alloc " +
                    $"{sizeGB} GB {sizeMB} MB {sizeKB} KB {sizeB} Byte.");
            }
            int size = (int)length;
            try
            {
                double[] data = GC.AllocateUninitializedArray<double>(size);
                return new NDArray(data, shape);
            }
            catch(Exception e)
            {
                throw new ArgumentOutOfRangeException(
                    "The array is too large! It need to alloc " +
                    $"{sizeGB} GB {sizeMB} MB {sizeKB} KB {sizeB} Byte.",
                    innerException: e);
            }
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

        public void Print(nint start = 6, nint end = 4,
            bool printMetadata = true, string? format = null)
        {
            using var _ = StringBuilderPool.Borrow(out var sb);
            ToString(sb, start, end, printMetadata, format);
            Console.WriteLine(sb.ToString());
        }

        public string ToString(nint start = 6, nint end = 4,
            bool printMetadata = true, string? format = null)
        {
            using var _ = StringBuilderPool.Borrow(out var sb);
            ToString(sb, start, end, printMetadata, format);
            return sb.ToString();
        }

        public void ToString(StringBuilder sb, nint start = 6,nint end = 4,
            bool printMetadata = true, string? format = null)
        {
            if(printMetadata)
                MetaDataString(sb);
            var segements = AsSegements();
            int dimLength = Math.Max(1, Rank);
            sb.Append('[', dimLength);
            if (segements.MoveNext())
            {
                segements.Current.ToString(sb, start, end, format);
                while (segements.MoveNext())
                {
                    sb.Append(']', segements.Step)
                    .Append(',')
                    .AppendLine()
                    .Append('[', segements.Step);
                    segements.Current.ToString(sb, start, end, format);
                }
            }
            sb.Append(']', dimLength);
        }

        public NDArraySegements AsSegements()
            => new(this);

        public NDArrayEnumerator GetEnumerator()
            => new(this);
        #endregion

        public ref struct NDArraySegements
        {
            private NDArray Array { get; }
            public nint Index { get; private set; }
            public nint Length { get; }
            public nint Batch { get; }
            public nint Stride { get; }
            public readonly Span<nint> StateSpan => state;

            public Span<nint> state;

            public ReadOnlySpan<nint> dimLengths;

            public ReadOnlySpan<nint> dimStrides;

            public VectorSpan Current { get; set; }

            public int Step { get; set; }

            public NDArraySegements(NDArray array)
            {
                if (array.Rank == 0)
                {
                    // For empty array, we can set it to the invalid state.
                    Array = array;
                    ReadOnlySpan<nint> lengths = array.Lengths;
                    Batch = lengths[0];
                    Length = 0;
                    Stride = 1;
                    state = [];
                    dimLengths = [];
                    dimStrides = [];
                    Index = 0;
                    Current = VectorSpan.Empty;
                }
                else if (array.Rank == 1)
                {
                    // For 1D array, we can use the VectorSpan directly.
                    Array = array;
                    ReadOnlySpan<nint> lengths = array.Lengths;
                    Batch = lengths[0];
                    Length = 1;
                    Stride = array.Strides[0];
                    state = [];
                    dimLengths = [];
                    dimStrides = [];
                    Index = -1;
                }
                else
                {
                    Array = array;
                    ReadOnlySpan<nint> lengths = array.Lengths;
                    Batch = lengths[^1];
                    Length = lengths[..^1].Product();
                    Stride = array.Strides[^1];
                    state = new nint[Array.Rank - 1];
                    var rank = Array.Rank;
                    dimLengths = Array.Metadata.AsSpan(0, rank - 1);
                    dimStrides = Array.Metadata.AsSpan(rank, rank - 1);
                    Index = -1;
                }
            }

            public bool MoveNext()
            {
                if(Length <= 1)
                {
                    if(Index >= 0)
                        return false;
                    else
                    {
                        Index = 0;
                        Current = new(ref Array.Data[Array.Offset], 
                            Batch, Stride);
                        return true;
                    }
                }
                int dimLength = dimLengths.Length;
                var stateSpan = StateSpan;
                if (Index > -1)
                {
                    Step = NintUtils.IncrementIndexLeft(dimLength - 1, stateSpan, dimLengths);
                }
                Index++;
                var stridesSpan = Array.Strides[..^1];
                nint index = Array.Offset + NintUtils.Dot(stateSpan, stridesSpan);
                ref double headRef = ref Array.Data[index];
                Current = new(ref headRef, Batch, Stride);
                return Index < Length;
            }
        }

        public ref struct NDArrayEnumerator(NDArray array)
        {
            NDArraySegements sequences = new(array);

            nint elementIndex = -1;

            public readonly ref double Current =>
                ref sequences.Current[elementIndex];

            public bool MoveNext()
            {
                elementIndex = (elementIndex + 1) % sequences.Batch;
                if (elementIndex == 0)
                    return sequences.MoveNext();
                return true;
            }
        }
    }
}
