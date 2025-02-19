using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace QuadMathSharp
{
    [StructLayout(LayoutKind.Sequential, Pack = 16, Size = 32)]
    public struct Complex128
    {
        public Float128 real;
        public Float128 imag;
    }
}
