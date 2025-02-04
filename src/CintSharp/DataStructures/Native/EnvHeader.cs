using NCDK.Numerics;
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
    [StructLayout(LayoutKind.Sequential, Size = 20)]
    public struct EnvHeader
    {

        /// <summary>
        /// Overall cutoff for integral prescreening, value needs to be ~ln(threshold)
        /// </summary>
        public double ExpCutoff { get; set; }

        /// <summary>
        /// R_C of (r-R_C) in dipole, GIAO operators
        /// </summary>
        public Vector3 CommonOrigin { get; set; }

        /// <summary>
        /// R_O in 1/|r-R_O|
        /// </summary>
        public Vector3 RinvOrigin { get; set; }

        /// <summary>
        /// ZETA parameter for Gaussian charge distribution (Gaussian nuclear model)
        /// </summary>
        public double RinvZeta { get; set; }

        /// <summary>
        /// omega parameter in range-separated coulomb operator
        /// LR interaction: erf(omega*r12)/r12 if omega > 0
        /// SR interaction: erfc(omega*r12)/r12 if omega < 0
        /// </summary>
        public double RangeOmega { get; set; }
        /// <summary>
        /// Yukawa potential and Slater-type geminal e^{-zeta r}
        /// </summary>
        public double F12Zeta { get; set; }
        /// <summary>
        /// Gaussian type geminal e^{-zeta r^2}
        /// </summary>
        public double GtgZeta { get; set; }
        /// <summary>
        /// Number of grids
        /// </summary>
        public double NumGrids { get; set; }
        /// <summary>
        /// 
        /// </summary>  
        public double Grids { get; set; }

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
