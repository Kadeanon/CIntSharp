﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace QuadMathSharp
{
    internal static class SoftFpExtentions
    {
        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static Float128 __addtf3(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static Float128 __subtf3(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static Float128 __multf3(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static Float128 __divtf3(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static Float128 __negtf2(Float128 a);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static int __eqtf2(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static int __netf2(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static int __lttf2(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static int __letf2(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static int __gttf2(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static int __getf2(Float128 a, Float128 b);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static Float128 __extendsftf2(float a);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static Float128 __extenddftf2(double a);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static float __trunctfsf2(Float128 a);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static double __trunctfdf2(Float128 a);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static double erf(double x);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static float erff(float x);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static double erfc(double x);

        [DllImport("libgcc_s_seh-1", CallingConvention = CallingConvention.Cdecl)]
        public extern static double erfcf(float x);
    }
}
