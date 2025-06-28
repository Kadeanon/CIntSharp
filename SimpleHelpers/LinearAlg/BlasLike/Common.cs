using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Text;
using System.Threading.Tasks;

namespace SimpleHelpers.LinearAlg
{
    public static partial class BlasLike
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Sum(this Vector256<double> vec)
        {
            return vec[0] + vec[1] + vec[2] + vec[3];
        }
    }
}
