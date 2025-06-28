using SimpleHelpers.MultiAlg.Helpers;
using SimpleHelpers.Utilities;
using System.Buffers;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {
        public ref double this[params ReadOnlySpan<nint> indexes]
        {
            get
            {
                int rank = Rank;
                int numIndex = indexes.Length;
                ArgumentOutOfRangeException.ThrowIfNotEqual(numIndex, rank,
                    $"The number of indexs ({numIndex}) must be " +
                    $"equal to the rank of the array ({rank}).");
                Span<NRange> newRanges = stackalloc NRange[rank];
                nint offset = Offset;
                int i = 0;
                for (; i < rank; i++)
                {
                    var dimOffset = indexes[i];
                    var dimTotalLength = Lengths[i];
                    ThrowUtils.ThrowIfNotInRange_LeftClosedRightOpen
                        (dimOffset, 0, dimTotalLength,
                        $"The index {dimOffset} is out of bounds " +
                        $"for dimension {i} with length {dimTotalLength}.");
                    offset += dimOffset * Strides[i];
                }
                return ref Data[offset];
            }
        }

        public ref double this[params ReadOnlySpan<NIndex> index]
        {
            get
            {
                int rank = Rank;
                int numIndex = index.Length;
                ArgumentOutOfRangeException.ThrowIfNotEqual(numIndex, rank,
                    $"The number of indexs ({numIndex}) must be " +
                    $"equal to the rank of the array ({rank}).");
                Span<NRange> newRanges = stackalloc NRange[rank];
                nint offset = Offset;
                int i = 0;
                for (; i < rank; i++)
                {
                    var dimIndex = index[i];
                    var dimTotalLength = Lengths[i];
                    var dimOffset =
                        dimIndex.GetOffset(dimTotalLength);
                    ThrowUtils.ThrowIfNotInRange_LeftClosedRightOpen
                        (dimOffset, 0, dimTotalLength,
                        $"The index {dimOffset} is out of bounds " +
                        $"for dimension {i} with length {dimTotalLength}.");
                    offset += dimOffset * Strides[i];
                }
                return ref Data[offset];
            }
        }

        public NDArray this[params ReadOnlySpan<NRange> ranges]
        {
            get => GetData(ranges);
            set => SetData(value, ranges);
        }

        public NDArray this[params ReadOnlySpan<NSlice> slices]
        {
            get => GetData(slices);
            set => SetData(value, slices);
        }

        public NDArray GetData(params ReadOnlySpan<NRange> ranges)
        {
            int rank = Rank;
            int numRanges = ranges.Length;
            ArgumentOutOfRangeException.ThrowIfGreaterThan(numRanges, rank,
                $"The number of ranges ({numRanges}) must be " +
                $"less than or equal to the rank of the array ({rank}).");
            nint offset = Offset;
            nint[] metadata = new nint[rank * 2];
            int i = 0;
            for (; i < ranges.Length; i++)
            {
                var dimRange = ranges[i];
                var dimTotalLength = Lengths[i];
                var (dimStart, dimLength) = ThrowUtils.ThrowIfNotInRange_CheckNRange
                    (dimRange, dimTotalLength, nameof(ranges));
                offset += dimStart * Strides[i];
                metadata[i] = dimLength;
                metadata[i + rank] = Strides[i];
            }
            for (; i < this.Rank; i++)
            {
                metadata[i] = Lengths[i];
                metadata[i + rank] = Strides[i];
            }
            return new NDArray(Data, offset, metadata, rank);
        }

        public void SetData(NDArray value, params ReadOnlySpan<NRange> ranges)
        {
            var array = GetData(ranges);
            SimpleTranspose.TransMethods.CopyTo(value, array);
        }

        public void SetData(double value, params ReadOnlySpan<NRange> ranges)
        {
            var array = GetData(ranges);
            array.Fill(value);
        }

        public NDArray GetData(params ReadOnlySpan<NSlice> slices)
        {
            int rank = Rank;
            var newRank = rank;
            int numSlices = slices.Length;
            ArgumentOutOfRangeException.ThrowIfGreaterThan(numSlices, rank,
                $"The number of ranges ({numSlices}) must be " +
                $"less than or equal to the rank of the array ({rank}).");
            Span<bool> indexers = stackalloc bool[slices.Length];
            Span<NRange> newRanges = stackalloc NRange[rank];
            nint offset = Offset;
            Span<nint> metadata = stackalloc nint[rank * 2];
            int i = 0;
            int j = 0;
            for (; i < numSlices; i++)
            {
                NSlice slice = slices[i];
                if (slice.IsIndex)
                {
                    var dimOffset = ThrowUtils.ThrowIfNotInRange_CheckNIndex
                        (slice.Index, Lengths[i], nameof(slice));
                    offset += dimOffset * Strides[i];
                }
                else
                {
                    var (dimOffset, dimLength) = ThrowUtils.ThrowIfNotInRange_CheckNRange
                        (slice.Range, Lengths[i], nameof(slice));
                    offset += dimOffset * Strides[i];
                    metadata[j] = dimLength;
                    metadata[j + rank] = Strides[i];
                    j++;
                }
            }
            if(j < i)
            {
                newRank -= i - j;
                for (int j2 = 0; j2 < j; j2++)
                {
                    metadata[j2 + newRank] = metadata[j2 + rank];
                }
                rank = newRank;
            }
            for (; i < rank; i++)
            {
                metadata[j] = Lengths[i];
                metadata[j + rank] = Strides[i + rank];
                j++;
            }
            return new NDArray(Data, offset, [ ..metadata[..(newRank * 2)]], newRank);
        }

        public void SetData(NDArray value, params ReadOnlySpan<NSlice> slices)
        {
            var array = GetData(slices);
            SimpleTranspose.TransMethods.CopyTo(value, array);
        }
    }
}
