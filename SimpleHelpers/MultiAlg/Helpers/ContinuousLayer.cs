using SimpleHelpers.Indices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace SimpleHelpers.MultiAlg.Helpers
{
    public unsafe struct ContinuousLayer
    {
        public bool IsHead { get; internal set; }
        public int Index { get; internal set; }
        public nint Start { get; internal set; }
        public nint Length { get; internal set; }
        public nint Size { get; internal set; }
        public nint Stride { get; internal set; }

        public static ContinuousLayer FromDiminfo(SingleIndice dim, int index)
        {
            return new ContinuousLayer()
            {
                IsHead = false,
                Index = index,
                Length = dim.Length,
                Size = dim.Length * dim.Stride,
                Stride = dim.Stride,
                Start = 0
            };
        }

        public static ContinuousLayer HeadFromDiminfo(SingleIndice dim, int index)
        {
            return new ContinuousLayer()
            {
                IsHead = true,
                Index = index,
                Length = dim.Length,
                Size = dim.Length * dim.Stride,
                Stride = dim.Stride,
                Start = 0
            };
        }

        internal bool TryAddDim(SingleIndice dim)
        {
            if (dim.Stride == Size)
            {
                Length *= dim.Length;
                Size *= dim.Length;
                return true;
            }
            return false;
        }
    }
}
