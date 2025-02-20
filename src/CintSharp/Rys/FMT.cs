using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using static CintSharp.Rys.Common;

namespace CintSharp.Rys
{
    internal unsafe static class FMT
    {
        /*
 * This code computes incomplete gamma function.  It is based on Xin
 * Wu's implementation.
 *
 *
 * List of Abbreviation(s)
 *
 * THO:
 * Gaussian-Expansion Methods for Molecular Integrals,
 * Hiroshi Taketa, Sigeru Huzinaga, and Kiyosi O-ohata,
 * Journal of the Physical Society of Japan,
 * Vol. 21, No. 11, 1966, 2313 - 2324.
 *
 */


        /*
         * Relative errors of fmt1_erfc_like are of
         *      (2*t)**(m-1) / (2m-3)!! * machine_precision * fmt_val
         * Errors of the other choice are
         *      (2m-1)!! / (2*t)**(m-1) * machine_precision * fmt_val
         * Given m, the turn-over point for t should satisfy
         *      (2m-1)!! / (2*t)**(m-1) > (2m-1)**.5
         * t0 = .5 * ((2m-1)!!/(2m-1)**.5)**(1/(m-1))
         */
        internal static readonly double[] TURNOVER_POINT = {
        0.0,
        0.0,
        0.866025403784,
        1.295010032056,
        1.705493613097,
        2.106432965305,
        2.501471934009,
        2.892473348218,
        3.280525047072,
        3.666320693281,
        4.05033123037 ,
        4.432891808508,
        4.814249856864,
        5.194593501454,
        5.574069276051,
        5.952793645111,
        6.330860773135,
        6.708347923415,
        7.08531930745 ,
        7.461828891625,
        7.837922483937,
        8.213639312398,
        8.589013237349,
        8.964073695432,
        9.338846443746,
        9.713354153046,
        10.08761688545,
        10.46165248270,
        10.83547688448,
        11.20910439128,
        11.58254788331,
        11.95581900374,
        12.32892831326,
        12.70188542111,
        13.07469909673,
        13.44737736550,
        13.81992759110,
        14.19235654675,
        14.56467047710,
        14.93687515212
};

        /*
         * Name
         *
         * fmtpse
         *
         * Synopsis
         *
         * double fmtpse(int m, double t)
         *
         * Description
         *
         * This function evaluates the auxiliary integral,
         *
         *             _ 1           2
         *            /     2 m  -t u
         * F (t)  =   |    u    e      du,
         *  m        _/  0
         *
         * by a power series expansion
         *
         *                    _                    2                     3                _
         *           Math.Exp(-t) |       t            t                     t                  |
         * F (t)  =  ------- | 1 + ----- + --------------- + ----------------------- + ... |
         *  m          2 b   |_    b + 1   (b + 1) (b + 2)   (b + 1) (b + 2) (b + 3)      _|,
         *
         * where b = m + 1 / 2. This power series expansion converges fast, when t is less than b + 1,
         * namely t < m + 3 / 2.
         *
         * Argument(s)
         *
         * int m:
         * F_m(t), see the Description section.
         *
         * double t:
         * F_m(t), see the Description section.
         *
         * Return Value
         * double:
         * F_m(t), see the Description section.
         *
         */

        /*
         * Name
         *
         * fmt
         *
         * Synopsis
         *
         * double fmt(int m, double t)
         *
         * Description
         *
         * This function evaluates the auxiliary integral, see Eq. 2.11 in THO,
         *
         *             _ 1           2
         *            /     2 m  -t u
         * F (t)  =   |    u    e      du,
         *  m        _/  0
         *
         * where m replaces ν for more convenient typesetting.
         *
         * If t is less than SML16 or equals 0, then
         *
         *              1
         * F (t)  =  -------.
         *  m        2 m + 1
         *
         * If t is less than m + 3 / 2, the auxiliary integral is evaluated by
         * a power series expansion (see fmtpse.c for details).
         *
         * Otherwise F (t) is calculated first
         *            0
         *                    1
         *                    -
         *           1 /  π  \2       _
         * F (t)  =  - | --- |  erf( /t ).
         *  0        2 \  t  /
         *
         * Then the upward recurrence relation is used for F (t) of higher m
         *                                                  m
         *
         *            (2 m - 1) F     (t) - Math.Exp( -t )
         *                       m - 1
         *  F (t)  =  -------------------------------.
         *   m                      2 t
         *
         * Argument(s)
         *
         * int m:
         * F_m(t), see the Description section.
         *
         * double t:
         * F_m(t), see the Description section.
         *
         * Return Value
         *
         * double:
         * F_m(t), see the Description section.
         *
         */
        internal static void fmt1_gamma_inc_like(double* f, double t, int m)
        {
            int i;
            double b = m + 0.5;
            double bi;
            double e = .5 * Math.Exp(-t);
            double x = e;
            double s = e;
            double tol = SML_FLOAT64 * e;
            for (bi = b + 1.0; x > tol; bi += 1.0)
            {
                x *= t / bi;
                s += x;
            }
            f[m] = s / b;
            for (i = m; i > 0; i--)
            {
                b -= 1.0;
                f[i - 1] = (e + t * f[i]) / b;
            }
        }

        internal static void gamma_inc_like(double* f, double t, int m)
        {
            if (t < TURNOVER_POINT[m])
            {
                fmt1_gamma_inc_like(f, t, m);
            }
            else
            {
                int i;
                double tt = Math.Sqrt(t);
                f[0] = (double)(SQRTPIE4 / tt * quad.Erf(tt));
                if (m > 0)
                {
                    double e = Math.Exp(-t);
                    double b = .5 / t;
                    for (i = 1; i <= m; i++)
                        f[i] = b * ((2 * i - 1) * f[i - 1] - e);
                }
            }
        }

        internal static void fmt1_lgamma_inc_like(quad* f, quad t, int m)
        {
            quad b = m + 0.5;
            quad bi;
            quad e = .5 * quad.Exp(-t);
            quad x = e;
            quad s = e;
            quad tol = SML_FLOAT80 * e;
            int i;
            for (bi = b + 1.0; x > tol; bi += 1.0)
            {
                x *= t / bi;
                s += x;
            }
            f[m] = s / b;
            for (i = m; i > 0; i--)
            {
                b -= 1;
                f[i - 1] = (e + t * f[i]) / b;
            }
        }

        internal static void lgamma_inc_like(quad* f, quad t, int m)
        {
            if (t < TURNOVER_POINT[m])
            {
                fmt1_lgamma_inc_like(f, t, m);
            }
            else
            {
                int i;
                quad tt = quad.Sqrt(t);
                f[0] = SQRTPIE4l / tt * quad.Erf(tt);
                if (m > 0)
                {
                    quad e = quad.Exp(-t);
                    quad b = .5 / t;
                    for (i = 1; i <= m; i++)
                        f[i] = b * ((2 * i - 1) * f[i - 1] - e);
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static double _pow(double @base, int exponent)
        {
            int i;
            double result = 1;
            for (i = 1; i <= exponent; i <<= 1)
            {
                if ((i & exponent) != 0)
                {
                    result *= @base;
                }
                @base *= @base;
            }
            return result;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static quad _powl(quad @base, int exponent)
        {
            int i;
            quad result = 1.0;
            for (i = 1; i <= exponent; i <<= 1)
            {
                if ((i & exponent) != 0)
                {
                    result *= @base;
                }
                @base *= @base;
            }
            return result;
        }

        /* This function evaluates the auxiliary integral,
         *
         *     2  _ 1           2
         *  t s  /     2 m  -t u
         * e     |    u    e      du,
         *      _/  s
         *
         * by a power series expansion
         *
         * F[m] = e^{t s^2} int_l^1 u^{2m} e^{-t u^2} du
         *      = e^{t s^2} /(2m+1) int e^{-t u^2} d u^{2m+1}
         *      = e^{t s^2} /(2m+1) [e^{-t u^2} u^{2m+1}]_l^1 + (2t)/(2m+1) int u^{2m+2} e^{-t u^2} du
         *      = e^{t s^2} /(m+.5) (.5*e^{-t} - .5*e^{-t l^2} l^{2m+1}) + t F[m+1])
         */
        internal static void fmt1_erfc_like(double* f, double t, double lower, int m)
        {
            int i;
            double lower2 = lower * lower;
            double b = m + 0.5;
            double bi;
            double e = .5 * Math.Exp(-t);
            double e1 = .5 * Math.Exp(-t * lower2) * lower;
            e1 *= _pow(lower2, m);
            double x = e;
            double x1 = e1;
            double s = e - e1;
            double div = 1.0;
            double delta = s;
            double tol = SML_FLOAT64 * Math.Abs(delta);
            for (bi = b + 1.0; Math.Abs(delta) > tol; bi += 1.0)
            {
                div *= t / bi;
                x1 *= lower2;
                delta = (x - x1) * div;
                s += delta;
            }
            double val = s / b;
            f[m] = val;
            for (i = m; i > 0; i--)
            {
                b -= 1.0;
                e1 /= lower2;
                val = (e - e1 + t * val) / b;
                f[i - 1] = val;
            }
        }
        internal static void fmt_erfc_like(double* f, double t, double lower, int m)
        {
            if (lower == 0)
            {
                gamma_inc_like(f, t, m);
                return;
            }

            int i;
            double lower2 = lower * lower;
            // F[m] < .5*sqrt(pi/t) * erfc(low*tt)
            if (t * lower2 > ERFC_bound)
            {
                for (i = 0; i <= m; i++)
                {
                    f[i] = 0;
                }
                return;
            }

            if (t < TURNOVER_POINT[m])
            {
                fmt1_erfc_like(f, t, lower, m);
            }
            else
            {
                double tt = Math.Sqrt(t);
                // erfc(a) - erfc(b) is more accurate than erf(b) - erf(a)
                double val = (double)(SQRTPIE4 / tt * (QuadMathSharp.SoftFpExtentions.erfc(lower * tt) - QuadMathSharp.SoftFpExtentions.erfc(tt)));
                f[0] = val;
                if (m > 0)
                {
                    double e = Math.Exp(-t);
                    double e1 = Math.Exp(-t * lower2) * lower;
                    double b = .5 / t;
                    for (i = 0; i < m; i++)
                    {
                        val = b * ((2 * i + 1) * val - e + e1);
                        e1 *= lower2;
                        f[i + 1] = val;
                    }
                }
            }
        }

        internal static void fmt_lerfc_like(quad* f, quad t, quad lower, int m)
        {
            if (lower == 0)
            {
                lgamma_inc_like(f, t, m);
                return;
            }

            int i;
            quad lower2 = lower * lower;
            // F[m] < .5*Math.Sqrt(pi/t) * erfc(low*tt)
            if (t * lower2 > ERFC_bound)
            {
                for (i = 0; i <= m; i++)
                {
                    f[i] = 0;
                }
                return;
            }

            if (t < TURNOVER_POINT[m])
            {
                fmt1_lerfc_like(f, t, lower, m);
            }
            else
            {
                quad tt = quad.Sqrt(t);
                // erfc(a) - erfc(b) is more accurate than erf(b) - erf(a)
                quad val = SQRTPIE4l / tt * (quad.Erfc(lower * tt) - quad.Erfc(tt));
                f[0] = val;
                if (m > 0)
                {
                    quad e = quad.Exp(-t);
                    quad e1 = quad.Exp(-t * lower2) * lower;
                    quad b = .5 / t;
                    for (i = 0; i < m; i++)
                    {
                        val = b * ((2 * i + 1) * val - e + e1);
                        e1 *= lower2;
                        f[i + 1] = val;
                    }
                }
            }
        }

        internal static void fmt1_lerfc_like(quad* f, quad t, quad lower, int m)
        {
            int i;
            quad lower2 = lower * lower;
            quad b = m + 0.5;
            quad bi;
            quad e = .5* quad.Exp(-t);
            quad e1 = .5* quad.Exp(-t * lower2) *lower;
            e1 *= _powl(lower2, m);
            quad x = e;
            quad x1 = e1;
            quad s = e - e1;
            quad div = 1.0;
            quad delta = s;
            quad tol = SML_FLOAT80 * quad.Abs(delta);
            for (bi = b + 1.0; quad.Abs(delta) > tol; bi += 1.0)
            {
                div *= t / bi;
                x1 *= lower2;
                delta = (x - x1) * div;
                s += delta;
            }
            quad val = s / b;
            f[m] = val;
            for (i = m; i > 0; i--)
            {
                b -= 1.0;
                e1 /= lower2;
                val = (e - e1 + t * val) / b;
                f[i - 1] = val;
            }
        }

        internal static void fmt1_qgamma_inc_like(__float128* f, __float128 t, int m)
        {
            __float128 b = m + .5;
            __float128 bi;
            __float128 e = .5 * quad.Exp(-t);
            __float128 x = e;
            __float128 s = e;
            __float128 tol = SML_FLOAT128 * e;
            int i;
            for (bi = b + 1.0; x > tol; bi += 1.0)
            {
                x *= t / bi;
                s += x;
            }
            f[m] = s / b;
            for (i = m; i > 0; i--)
            {
                b -= 1;
                f[i - 1] = (e + t * f[i]) / b;
            }
        }

        internal static void qgamma_inc_like(__float128* f, __float128 t, int m)
        {
            if (t < TURNOVER_POINT[m])
            {
                fmt1_qgamma_inc_like(f, t, m);
            }
            else
            {
                int i;
                __float128 tt = quad.Sqrt(t);
                f[0] = SQRTPIE4q / tt * quad.Erf(tt);
                if (m > 0)
                {
                    __float128 e = quad.Exp(-t);
                    __float128 b = 0.5 / t;
                    for (i = 1; i <= m; i++)
                        f[i] = b * ((2 * i - 1) * f[i - 1] - e);
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static  __float128 _powq(__float128 @base, int exponent)
        {
            int i;
            __float128 result = 1.1;
            for (i = 1; i <= exponent; i <<= 1)
            {
                if ((i & exponent) != 0)
                {
                    result *= @base;
                }
                @base *= @base;
            }
            return result;
        }

        internal static void fmt_qerfc_like(__float128* f, __float128 t, __float128 lower, int m)
        {
            if (lower == 0)
            {
                qgamma_inc_like(f, t, m);
                return;
            }

            int i;
            __float128 lower2 = lower * lower;
            // F[m] < .5*Math.Sqrt(pi/t) * erfc(low*tt)
            if (t * lower2 > ERFC_bound)
            {
                for (i = 0; i <= m; i++)
                {
                    f[i] = 0;
                }
                return;
            }

            if (t < TURNOVER_POINT[m])
            {
                fmt1_qerfc_like(f, t, lower, m);
            }
            else
            {
                __float128 tt = quad.Sqrt(t);
                // erfc(a) - erfc(b) is more accurate than erf(b) - erf(a)
                __float128 val = SQRTPIE4q / tt * (quad.Erfc(lower * tt) - quad.Erfc(tt));
                f[0] = val;
                if (m > 0)
                {
                    __float128 e = quad.Exp(-t);
                    __float128 e1 = quad.Exp(-t * lower2) * lower;
                    __float128 b = 0.5 / t;
                    for (i = 0; i < m; i++)
                    {
                        val = b * ((2 * i + 1) * val - e + e1);
                        e1 *= lower2;
                        f[i + 1] = val;
                    }
                }
            }
        }

        internal static void fmt1_qerfc_like(__float128* f, __float128 t, __float128 lower, int m)
        {
            int i;
            __float128 lower2 = lower * lower;
            __float128 b = m + .5;
            __float128 bi;
            __float128 e = .5 * quad.Exp(-t);
            __float128 e1 = .5 * quad.Exp(-t * lower2) *lower;
            e1 *= _powq(lower2, m);
            __float128 x = e;
            __float128 x1 = e1;
            __float128 s = e - e1;
            __float128 div = 1.0;
            __float128 delta = s;
            __float128 tol = SML_FLOAT128 * quad.Abs(delta);
            for (bi = b + 1.0; quad.Abs(delta) > tol; bi += 1.0)
            {
                div *= t / bi;
                x1 *= lower2;
                delta = (x - x1) * div;
                s += delta;
            }
            __float128 val = s / b;
            f[m] = val;
            for (i = m; i > 0; i--)
            {
                b -= 1.0;
                e1 /= lower2;
                val = (e - e1 + t * val) / b;
                f[i - 1] = val;
            }
        }

    }
}
