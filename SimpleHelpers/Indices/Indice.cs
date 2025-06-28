global using LengthType = nint;
global using LengthIndex = System.Buffers.NIndex;
global using NRange = System.Buffers.NRange;
global using StrideType = nint;
global using OffsetType = nint;

namespace SimpleHelpers.Indices
{
    public readonly record struct SingleIndice(LengthType Length, StrideType Stride)
    { }

    public readonly record struct DoubleIndice(LengthType Length, StrideType AStride,
        StrideType BStride)
    {
        public SingleIndice A => new(Length, AStride);

        public SingleIndice B => new(Length, BStride);

        public DoubleIndice Swap() => new(Length, BStride, AStride);
    }

    public readonly record struct TripleIndice(LengthType Length, StrideType AStride,
        StrideType BStride, StrideType CStride)
    {
        public SingleIndice A => new(Length, AStride);
        public SingleIndice B => new(Length, BStride);
        public SingleIndice C => new(Length, CStride);

        public DoubleIndice AB => new(Length, AStride, BStride);
        public DoubleIndice AC => new(Length, AStride, CStride);
        public DoubleIndice BC => new(Length, BStride, CStride);
    }
}
