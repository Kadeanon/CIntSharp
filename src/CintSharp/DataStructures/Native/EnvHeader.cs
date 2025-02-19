using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Net.NetworkInformation;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DataStructures.Native
{
    [StructLayout(LayoutKind.Explicit, Size = 20)]
    public struct EnvHeader
    {

        public struct Buffer
        {
            public double ofset;
        }

        /// <summary>
        /// Overall cutoff for integral prescreening, value needs to be ~ln(threshold)
        /// </summary>
        public ref double ExpCutoff => ref Span[0];

        /// <summary>
        /// R_C of (r-R_C) in dipole, GIAO operators
        /// </summary>
        public ref Vector3 CommonOrigin => ref Unsafe.As<double, Vector3>(ref Span[1]);

        /// <summary>
        /// R_O in 1/|r-R_O|
        /// </summary>
        public ref Vector3 RinvOrigin => ref Unsafe.As<double, Vector3>(ref Span[4]);

        /// <summary>
        /// ZETA parameter for Gaussian charge distribution (Gaussian nuclear model)
        /// </summary>
        public ref double RinvZeta => ref Span[7];

        /// <summary>
        /// omega parameter in range-separated coulomb operator
        /// LR interaction: erf(omega*r12)/r12 if omega > 0
        /// SR interaction: erfc(omega*r12)/r12 if omega < 0
        /// </summary>
        public ref double RangeOmega => ref Span[8];
        /// <summary>
        /// Yukawa potential and Slater-type geminal e^{-zeta r}
        /// </summary>
        public ref double F12Zeta => ref Span[9];
        /// <summary>
        /// Gaussian type geminal e^{-zeta r^2}
        /// </summary>
        public ref double GtgZeta => ref Span[10];
        /// <summary>
        /// Number of grids
        /// </summary>
        public ref double NumGrids => ref Span[11];
        /// <summary>
        /// 
        /// </summary>  
        public ref double Grids => ref Span[12];

        public const int Length = 13;
        public const int PointerEnvStart = 20;

        public void WriteTo(Span<double> arr) 
        {
            Span.CopyTo(arr);
        }

        public Span<double> Span => 
            MemoryMarshal.CreateSpan(ref Unsafe.As<EnvHeader, double>(ref this), Length);

        public static ref EnvHeader FromSpan(Span<double> span)
        {
            return ref Unsafe.As<double, EnvHeader>(ref MemoryMarshal.GetReference(span));
        }
    }
}
