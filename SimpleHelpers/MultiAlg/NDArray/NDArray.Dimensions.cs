using SimpleHelpers.Indices;
using SimpleHelpers.MultiAlg.Helpers;
using SimpleHelpers.Utilities;
using System.Diagnostics.CodeAnalysis;
using System.Numerics.Tensors;
using System.Runtime.InteropServices;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {

        public NDArray Flatten(bool forceCopy = false)
        {
            if (Rank == 1 && !forceCopy)
            {
                return this;
            }
            nint[] newShape = { Size };
            return Reshape(newShape, forceCopy);
        }

        public NDArray Reshape(params ReadOnlySpan<nint> shape) => Reshape(shape, true, false);

        public NDArray Reshape(ReadOnlySpan<nint> shape, bool allowCopy = true, bool forceCopy = false)
        {
            Span<nint> tempShape = stackalloc nint[shape.Length];
            shape.CopyTo(tempShape);
            //如果形状相同，直接返回
            if (Reshape_CheckShape(tempShape) && !forceCopy)
            {
                return this;
            }

            var view = Reshape_TryBuildHead(tempShape, out var head);
            if (!view || forceCopy)
            {
                if (!allowCopy)
                {
                    throw new InvalidOperationException("Cannot reshape without copying");
                }
                return Reshape_Copy(tempShape);
            }
            else
            {
                return new(Data, Offset, head, tempShape.Length);
            }
        }

        private bool Reshape_CheckShape(Span<nint> shape)
        {
            if (shape.Length == 0)
            {
                throw new ArgumentException("Shape cannot be empty");
            }
            nint size = 1;
            int wildcard = -1;
            for (int i = 0; i < shape.Length; i++)
            {
                nint shapeValue = shape[i];
                if (shapeValue <= 0)
                {
                    if (shapeValue != -1)
                    {
                        throw new ArgumentException("Shape values must be positive or -1 as wildcard");
                    }
                    if (wildcard != -1)
                    {
                        throw new ArgumentException("Shape cannot contain more than one wildcard");
                    }
                    wildcard = i;
                    continue;
                }
                size *= shape[i];
            }
            if (wildcard != -1)
            {
                (var quotient, var remainder) = Math.DivRem(Size, size);
                if (remainder != 0)
                {
                    throw new ArgumentException("Cannot determine wildcard value, the size of array is not a multiple of the new size with other dims.");
                }
                shape[wildcard] = quotient;
            }
            else
            {
                if (size != Size)
                {
                    throw new ArgumentException("Size of new shape must be equal to the size of the original shape");
                }
            }
            return shape.Length == Lengths.Length && shape.SequenceEqual(Lengths);
        }

        private bool Reshape_TryBuildHead(ReadOnlySpan<nint> shape, out nint[] head)
        {
            var rank = shape.Length;
            head = new nint[rank * 2];
            ContinuousInfo info = new(Metadata, false);
            int newNumDims = shape.Length;
            Span<nint> strides = stackalloc nint[newNumDims];
            Span<nint> starts = stackalloc nint[newNumDims];
            //Try check with continuous info
            Span<ContinuousLayer> Continuouslayers = CollectionsMarshal.AsSpan(info.Layers);
            Span<int> layers = stackalloc int[Continuouslayers.Length];
            int numLayers = 0;
            for (var iLayer = 0; iLayer < Continuouslayers.Length; iLayer++)
            {
                ref var layer = ref Continuouslayers[iLayer];
                if (layer.IsHead)
                {
                    layers[numLayers] = iLayer;
                    numLayers++;
                }
            }
            layers = layers[..numLayers];
            numLayers--;
            ref ContinuousLayer layerHeader = ref Continuouslayers[layers[numLayers]];
            nint stride = layerHeader.Stride;
            nint size = 1;
            nint start = layerHeader.Start;
            for (int i = newNumDims - 1; i >= 0; i--)
            {
                strides[i] = stride;
                starts[i] = start;
                stride *= shape[i];
                start = 0;
                size *= shape[i];
                if (size > layerHeader.Size)
                {
                    return false;
                }
                else if (size == layerHeader.Size)
                {
                    if (numLayers == 0)
                    {
                        break;
                    }
                    numLayers--;
                    layerHeader = ref Continuouslayers[layers[numLayers]];
                    stride = layerHeader.Stride;
                    size = 1;
                    start = layerHeader.Start;
                }
            }
            for (int i = 0; i < newNumDims; i++)
            {
                head[i] = shape[i];
                head[i + rank] = strides[i];
            }
            return true;
        }

        private NDArray Reshape_Copy(Span<nint> shape)
        {
            double[] arrayData = new double [TensorPrimitives.Product<nint>(shape)];
            FlattenTo(arrayData);
            return new(arrayData, shape);
        }

        public bool TryReshapeWithoutCopy(ReadOnlySpan<nint> shape, 
            [NotNullWhen(true)] out NDArray? result)
        {
            Span<nint> tempShape = stackalloc nint[shape.Length];
            shape.CopyTo(tempShape);
            if (Reshape_CheckShape(tempShape))
                result = this;

            if (Reshape_TryBuildHead(tempShape, out var head))
            {
                result = new(Data, Offset, head, tempShape.Length);
                return true;
            }
            else
            {
                result = null;
                return false;
            }
        }

        public NDArray SwapAxis(int dim0 = 0, int dim1 = 1)
        {
            var newArray = View();
            (Metadata[dim0], Metadata[dim1]) =
                (Metadata[dim1], Metadata[dim0]);
            (Metadata[Rank + dim0], Metadata[Rank + dim1]) =
                (Metadata[Rank + dim1], Metadata[Rank + dim0]);
            return newArray;
        }

        public NDArray Squeeze()
        {
            Span<SingleIndice> indices = stackalloc SingleIndice[Rank];
            int newRank = 0;
            var lengths = Lengths;
            var strides = Strides;
            for (int i = 0; i < indices.Length; i++)
            {
                var length = lengths[i];
                var stride = strides[i];
                if (indices[i].Length != 1)
                {
                    indices[newRank] = indices[i];
                    newRank++;
                }
            }
            if (newRank == Rank)
            {
                return this; // No change needed
            }
            else
            {
                nint[] newShape = new nint[newRank * 2];
                var newLengths = newShape.AsSpan(0, newRank);
                var newStrides = newShape.AsSpan(newRank, newRank);
                for (int i = 0; i < newRank; i++)
                {
                    newLengths[i] = indices[i].Length;
                    newStrides[i] = indices[i].Stride;
                }
                return new NDArray(Data, Offset, newShape, newRank);
            }
        }

        public NDArray SqueezeDimension(int dim)
        {
            if (dim < 0 || dim >= Rank)
                throw new ArgumentOutOfRangeException(nameof(dim),
                    $"The dim value {dim} is out of range [0, {Rank})");
            var length = Lengths[dim];
            if (length != 1)
                throw new ArgumentException(
                    $"The dimension {dim} has length {length}, " +
                    "which cannot be squeezed.", nameof(dim));
            int newRank = Rank - 1;
            nint[] newShape = new nint[newRank * 2];
            Array.Copy(Metadata, newShape, dim);
            //Skip Lengths[dim]
            Array.Copy(Metadata, dim + 1, newShape, dim, Rank - 1);
            //Skip Strides[dim]
            Array.Copy(Metadata, Rank + dim + 1, newShape,
                newRank + dim, Rank - dim - 1);
            return new NDArray(Data, Offset, newShape, newRank);
        }

        public NDArray SliceFirstDim(nint index)
        {
            nint[] metadata = new nint[(Rank - 1) * 2];
            nint length0 = Lengths[0];
            nint stride0 = Strides[0];
            ThrowUtils.ThrowIfNotInRange_CheckNegative
                (ref index, length0,
                $"The index {index} is out of bounds for the first dimension with length {length0}.");
            nint offset = Offset + index * stride0;
            for (int i = 1; i < Rank; i++)
            {
                metadata[i - 1] = Lengths[i];
                metadata[Rank + i - 2] = Strides[i];
            }
            return new NDArray(Data, offset, metadata, Rank - 1);
        }

        public NDArray Transpose(params ReadOnlySpan<int> dims)
        {
            var newArray = View();
            var head = newArray.Metadata;
            Span<bool> used = stackalloc bool[Rank];
            if (dims.Length == Rank)
            {
                var lengths = Lengths;
                var strides = Strides;
                var newLengths = newArray.Metadata.AsSpan(0, Rank);
                var newStrides = newArray.Metadata.AsSpan(Rank);
                for (int i = 0; i < Rank; i++)
                {
                    var dim = dims[i];
                    if (dim < 0 || dim >= Rank)
                        throw new ArgumentOutOfRangeException(nameof(dims),
                            $"The dim value {dim} is out of range [0, {Rank})");
                    if (used[dim])
                        throw new ArgumentException(nameof(dims),
                            "The axis must not contain duplicate values");
                    used[dim] = true;
                    newLengths[i] = lengths[dim];
                    newStrides[i] = strides[dim];
                }
                return newArray;
            }
            else if (dims.Length == 0)
            {
                var newLengths = newArray.Metadata.AsSpan(0, Rank);
                var newStrides = newArray.Metadata.AsSpan(Rank);
                newLengths.Reverse();
                newStrides.Reverse();
                return newArray;
            }
            throw new ArgumentException("The length of the axis must be equal to the rank of the array or empty", nameof(dims));
        }

        public NDArray View()
        {
            return new NDArray(Data,
                Offset, (nint[])Metadata.Clone(), Rank);
        }
    }
}
