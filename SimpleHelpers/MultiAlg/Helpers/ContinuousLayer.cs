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
        (SingleIndice dim, int index, bool isHead = false)
    {
        public bool IsHead { get; internal set; } = isHead;
        public int Index { get; internal set; } = index;
        public nint Length { get; internal set; } = dim.Length;
        public nint BlockSize { get; internal set; } = dim.Length * dim.Stride;
        public nint Stride { get; internal set; } = dim.Stride;

        public override readonly string ToString()
        {
            if(IsHead)
            {
                return $"Head Layer {Index}: Length={Length}, Size={BlockSize}, Stride={Stride}";
            }
            else
            {
                return $"Layer {Index}: Length={Length}, Size={BlockSize}, Stride={Stride}";
            }
        }
    }
}
