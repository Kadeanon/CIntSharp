using System.Buffers;

namespace SimpleHelpers.MultiAlg.Helpers
{
    public struct NSlice
    {
        NIndex index;
        NRange range;
        bool isIndex;

        public NSlice(NIndex index)
        {
            this.index = index;
            range = default;
            isIndex = true;
        }

        public static implicit operator NSlice(NIndex index)
            => new NSlice(index);

        public static implicit operator NSlice(NRange range)
            => new NSlice(range);

        public static implicit operator NSlice(Index index)
            => new NSlice(new NIndex(index));

        public static implicit operator NSlice(Range range)
            => new NSlice(new NRange(range));

        public static implicit operator NSlice(int index)
            => new NSlice(new NIndex(index));

        public static implicit operator NSlice(nint index)
            => new NSlice(new NIndex(index));

        public NIndex Index => index;

        public NRange Range => range;

        public bool IsIndex => isIndex;

        public NSlice(NRange range)
        {
            index = default;
            this.range = range;
            isIndex = false;
        }
    }
}
