global using quad = QuadMathSharp.Float128;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace CintSharp.Rys
{
    internal static class Common
    {
        internal const double PIE4 = 0.78539816339744827900;
        internal const double DBL_EPSILON = 2.2204460492503131E-16;
        internal const double THRESHOLD_ZERO = DBL_EPSILON * 8;
        internal const double SMALLX_LIMIT = 3e-7;

        // other boundaries
        internal const int ANG_MAX = 15; // l = 0..15
        internal const int LMAX1 = 16; // > ANG_MAX
        internal const int CART_MAX = 136; // > (ANG_MAX*(ANG_MAX+1)/2)
        internal const int SHLS_MAX = 1048576;
        internal const int NPRIM_MAX = 64;
        internal const int NCTR_MAX = 64;

        internal const int POINT_NUC = 1;
        internal const int GAUSSIAN_NUC = 2;
        internal const int FRAC_CHARGE_NUC = 3;



        internal static readonly quad SQRTPIE4 = quad.FromString("0.8862269254527580136490837416705725913987747280611935641069038949264");
        internal static readonly quad SQRTPIE4l = quad.FromString("0.8862269254527580136490837416705725913987747280611935641069038949264");


        internal static readonly quad SQRTPIE4q = quad.FromString("0.8862269254527580136490837416705725913987747280611935641069038949264");

        // For double precision 0.063^20 < 2**-53
        internal const int FLOCKE_EXTRA_ORDER_FOR_DP = 20;
        // For long double 0.063^24 < 2**-64
        internal const int FLOCKE_EXTRA_ORDER_FOR_LP = 24;
        // For quadruple precision 0.063^36 < 2**-113
        internal const int FLOCKE_EXTRA_ORDER_FOR_QP = 36;

        internal const double SML_FLOAT64 = (DBL_EPSILON * .5);
        internal const double SML_FLOAT80 = 2.0e-20;
        internal const int ERFC_bound = 200;
        internal const double SML_FLOAT128 = 1.0e-35;
    }
}
