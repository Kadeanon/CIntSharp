using CintSharp.DataStructures;
using CintSharp.DataStructures.Native;
using CintSharp.Rys;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Intors.Int1e
{
    internal unsafe static class G1e
    {
        internal static void CINTg1e_index_xyz(FINT[] idx, CINTEnvVars envs)
        {
            FINT i_l = envs.i_l;
            FINT j_l = envs.j_l;
            FINT nfi = envs.nfi;
            FINT nfj = envs.nfj;
            FINT di = envs.g_stride_i;
            FINT dj = envs.g_stride_j;
            FINT i, j, n;
            FINT ofx, ofjx;
            FINT ofy, ofjy;
            FINT ofz, ofjz;
            Span<FINT> i_nx = stackalloc FINT[CART_MAX],
                i_ny = stackalloc FINT[CART_MAX],
                i_nz = stackalloc FINT[CART_MAX],
                j_nx = stackalloc FINT[CART_MAX],
                j_ny = stackalloc FINT[CART_MAX],
                j_nz = stackalloc FINT[CART_MAX];

            CINTcart_comp(i_nx, i_ny, i_nz, i_l);
            CINTcart_comp(j_nx, j_ny, j_nz, j_l);

            ofx = 0;
            ofy = envs.g_size;
            ofz = envs.g_size * 2;
            n = 0;
            for (j = 0; j < nfj; j++)
            {
                ofjx = ofx + dj * j_nx[j];
                ofjy = ofy + dj * j_ny[j];
                ofjz = ofz + dj * j_nz[j];
                for (i = 0; i < nfi; i++)
                {
                    idx[n + 0] = ofjx + di * i_nx[i];
                    idx[n + 1] = ofjy + di * i_ny[i];
                    idx[n + 2] = ofjz + di * i_nz[i];
                    n += 3;
                }
            }
        }

        internal static FINT CINTg1e_ovlp(Span<double> g, CINTEnvVars envs)
        {
            Debug.Assert(g.Length >= envs.g_size * 3);
            FINT size = envs.g_size;
            Span<double> gx = g[..size];
            Span<double> gy = g[size..(size * 2)];
            Span<double> gz = g[(size * 2)..];
            double aij = envs.ai + envs.aj;

            gx[0] = 1;
            gy[0] = 1;
            gz[0] = envs.fac * Math.Sqrt(Math.PI) * Math.PI / (aij * Math.Sqrt(aij));

            FINT nmax = envs.li_ceil + envs.lj_ceil;
            if (nmax == 0)
            {
                return 1;
            }

            Vector3 rij = envs.rij;
            Vector3 rirj = envs.rirj;
            FINT lj, di, dj;
            FINT i, j, n, ptr;
            Vector3 rx;
            if (envs.li_ceil > envs.lj_ceil)
            {
                // li = envs.li_ceil;
                lj = envs.lj_ceil;
                di = envs.g_stride_i;
                dj = envs.g_stride_j;
                rx = envs.ri;
            }
            else
            {
                // li = envs.lj_ceil;
                lj = envs.li_ceil;
                di = envs.g_stride_j;
                dj = envs.g_stride_i;
                rx = envs.rj;
            }
            Vector3 rijrx = rij - rx;

            gx[di] = rijrx[0] * gx[0];
            gy[di] = rijrx[1] * gy[0];
            gz[di] = rijrx[2] * gz[0];

            double aij2 = .5 / aij;
            for (i = 1; i < nmax; i++)
            {
                gx[(i + 1) * di] = i * aij2 * gx[(i - 1) * di] + rijrx[0] * gx[i * di];
                gy[(i + 1) * di] = i * aij2 * gy[(i - 1) * di] + rijrx[1] * gy[i * di];
                gz[(i + 1) * di] = i * aij2 * gz[(i - 1) * di] + rijrx[2] * gz[i * di];
            }

            for (j = 1; j <= lj; j++)
            {
                ptr = dj * j;
                for (i = 0, n = ptr; i <= nmax - j; i++, n += di)
                {
                    gx[n] = gx[n + di - dj] + rirj[0] * gx[n - dj];
                    gy[n] = gy[n + di - dj] + rirj[1] * gy[n - dj];
                    gz[n] = gz[n + di - dj] + rirj[2] * gz[n - dj];
                }
            }
            return 1;
        }

        /*
         * Calculate temporary parameter tau for nuclear charge distribution.
         * The charge parameter zeta is defined as    rho(r) = Norm * exp(-zeta*r^2)
         */
        internal static double CINTnuc_mod(double aij, FINT nuc_id, Atm[] atm, double[] env)
        {
            double zeta;
            if (nuc_id < 0)
            {
                zeta = env[PTR_RINV_ZETA];
            }
            else if (atm[nuc_id].NuclearModel == NuclearModelKind.GaussianNuclear)
            {
                zeta = env[atm[nuc_id].pointerForZeta];
            }
            else
            {
                zeta = 0;
            }

            if (zeta > 0)
            {
                return Math.Sqrt(zeta / (aij + zeta));
            }
            else
            {
                return 1;
            }
        }

        internal static FINT CINTg1e_nuc(Span<double> g, CINTEnvVars envs, FINT nuc_id)
        {
            FINT nrys_roots = envs.nrys_roots;
            Atm[] atm = envs.atm;
            double[] env = envs.env;
            Vector3 rij = envs.rij;
            Debug.Assert(g.Length >= envs.g_size * 3);
            FINT size = envs.g_size;
            Span<double> gx = g[..size];
            Span<double> gy = g[size..(size * 2)];
            Span<double> gz = g[(size * 2)..];
            Span<double> u = stackalloc double[MXRYSROOTS];
            Span<double> w = gz;
            ref Vector3 cr = ref Unsafe.NullRef<Vector3>();
            FINT i, j, n;
            Vector3 crij;
            double x, fac1;
            double aij = envs.ai + envs.aj;
            double tau = CINTnuc_mod(aij, nuc_id, atm, env);
            ref EnvHeader header = ref EnvHeader.FromSpan(env);
            if (nuc_id < 0)
            {
                fac1 = 2 * Math.PI * envs.fac * tau / aij;
                cr = ref header.RinvOrigin;
            }
            else if (atm[nuc_id].NuclearModel == NuclearModelKind.FracChargeNuclear)
            {
                fac1 = 2 * Math.PI * - env[atm[nuc_id].pointerForFracCharge] * envs.fac * tau / aij;
                cr = ref Unsafe.As<double, Vector3>(ref env[atm[nuc_id].pointerForCoords]);
            }
            else
            {
                fac1 = 2 * Math.PI * -Math.Abs(atm[nuc_id].charge) * envs.fac * tau / aij;
                cr = ref Unsafe.As<double, Vector3>(ref env[atm[nuc_id].pointerForCoords]);
            }
            crij = cr - rij;
            x = aij * tau * tau * crij[0] *crij[0];
            RysRoot.CINTrys_roots(nrys_roots, x, u, w);

            for (i = 0; i < nrys_roots; i++)
            {
                gx[i] = 1;
                gy[i] = 1;
                gz[i] *= fac1;
            }
            FINT nmax = envs.li_ceil + envs.lj_ceil;
            if (nmax == 0)
            {
                return 1;
            }

            double*
                p0x,
                p0y,
                p0z,
                p1x,
                p1y,
                p1z,
                p2x,
                p2y,
                p2z;
            FINT lj, di, dj;
            ref Vector3 rx = ref Unsafe.NullRef<Vector3>();
            if (envs.li_ceil > envs.lj_ceil)
            {
                // li = envs.li_ceil;
                lj = envs.lj_ceil;
                di = envs.g_stride_i;
                dj = envs.g_stride_j;
                rx = ref envs.ri;
            }
            else
            {
                // li = envs.lj_ceil;
                lj = envs.li_ceil;
                di = envs.g_stride_j;
                dj = envs.g_stride_i;
                rx = ref envs.rj;
            }
            double rijrx = rij[0] - rx[0];
            double rijry = rij[1] - rx[1];
            double rijrz = rij[2] - rx[2];
            double aij2 = 0.5 / aij;
            double ru, rt, r0, r1, r2;

            p0x = (double*)Unsafe.AsPointer(ref gx[di]);
            p0y = (double*)Unsafe.AsPointer(ref gy[di]);
            p0z = (double*)Unsafe.AsPointer(ref gz[di]);
            p1x = (double*)Unsafe.AsPointer(ref gx[di]);
            p1y = (double*)Unsafe.AsPointer(ref gy[di]);
            p1z = (double*)Unsafe.AsPointer(ref gz[di]);
            for (n = 0; n < nrys_roots; n++)
            {
                ru = tau * tau * u[n] / (1 + u[n]);
                rt = aij2 - aij2 * ru;
                r0 = rijrx + ru * crij[0];
                r1 = rijry + ru * crij[1];
                r2 = rijrz + ru * crij[2];

                p0x[n] = r0 * gx[n];
                p0y[n] = r1 * gy[n];
                p0z[n] = r2 * gz[n];
                for (i = 1; i < nmax; i++)
                {
                    p0x[n + i * di] = i * rt * p1x[n + i * di] + r0 * gx[n + i * di];
                    p0y[n + i * di] = i * rt * p1y[n + i * di] + r1 * gy[n + i * di];
                    p0z[n + i * di] = i * rt * p1z[n + i * di] + r2 * gz[n + i * di];
                }
            }

            double rirjx = envs.rirj[0];
            double rirjy = envs.rirj[1];
            double rirjz = envs.rirj[2];
            for (j = 1; j <= lj; j++)
            {
                p0x = (double*)Unsafe.AsPointer(ref gx[j * dj]);
                p0y = (double*)Unsafe.AsPointer(ref gy[j * dj]);
                p0z = (double*)Unsafe.AsPointer(ref gz[j * dj]);
                p1x = p0x - dj;
                p1y = p0y - dj;
                p1z = p0z - dj;
                p2x = p1x + di;
                p2y = p1y + di;
                p2z = p1z + di;
                for (i = 0; i <= nmax - j; i++)
                {
                    for (n = 0; n < nrys_roots; n++)
                    {
                        p0x[n + i * di] = p2x[n + i * di] + rirjx * p1x[n + i * di];
                        p0y[n + i * di] = p2y[n + i * di] + rirjy * p1y[n + i * di];
                        p0z[n + i * di] = p2z[n + i * di] + rirjz * p1z[n + i * di];
                    }
                }
            }
            return 1;
        }

        internal static void CINTnabla1i_1e(double* f, double* g,
                            FINT li, FINT lj, FINT lk, CINTEnvVars envs)
        {
            FINT dj = envs.g_stride_j;
            FINT dk = envs.g_stride_k;
            double ai2 = -2 * envs.ai;
            FINT i, j, k, ptr;
            double* gx = g;
            double* gy = g + envs.g_size;
            double* gz = g + envs.g_size * 2;
            double* fx = f;
            double* fy = f + envs.g_size;
            double* fz = f + envs.g_size * 2;

            for (k = 0; k <= lk; k++)
            {
                for (j = 0; j <= lj; j++)
                {
                    ptr = dj * j + dk * k;
                    //f(...,0,...) = -2*ai*g(...,1,...)
                    fx[ptr] = ai2 * gx[ptr + 1];
                    fy[ptr] = ai2 * gy[ptr + 1];
                    fz[ptr] = ai2 * gz[ptr + 1];
                    //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                    for (i = 1; i <= li; i++)
                    {
                        fx[ptr + i] = i * gx[ptr + i - 1] + ai2 * gx[ptr + i + 1];
                        fy[ptr + i] = i * gy[ptr + i - 1] + ai2 * gy[ptr + i + 1];
                        fz[ptr + i] = i * gz[ptr + i - 1] + ai2 * gz[ptr + i + 1];
                    }
                }
            }
        }

        internal static void CINTnabla1j_1e(double* f, double* g,
                            FINT li, FINT lj, FINT lk, CINTEnvVars envs)
        {
            FINT dj = envs.g_stride_j;
            FINT dk = envs.g_stride_k;
            double aj2 = -2 * envs.aj;
            FINT i, j, k, ptr;
            double* gx = g;
            double* gy = g + envs.g_size;
            double* gz = g + envs.g_size * 2;
            double* fx = f;
            double* fy = f + envs.g_size;
            double* fz = f + envs.g_size * 2;

            for (k = 0; k <= lk; k++)
            {
                ptr = dk * k;
                //f(...,0,...) = -2*aj*g(...,1,...)
                for (i = ptr; i <= ptr + li; i++)
                {
                    fx[i] = aj2 * gx[i + dj];
                    fy[i] = aj2 * gy[i + dj];
                    fz[i] = aj2 * gz[i + dj];
                }
                //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
                for (j = 1; j <= lj; j++)
                {
                    ptr = dj * j + dk * k;
                    for (i = ptr; i <= ptr + li; i++)
                    {
                        fx[i] = j * gx[i - dj] + aj2 * gx[i + dj];
                        fy[i] = j * gy[i - dj] + aj2 * gy[i + dj];
                        fz[i] = j * gz[i - dj] + aj2 * gz[i + dj];
                    }
                }
            }
        }

        /*
         * ( ij | \nabla k )
         */
        internal static void CINTnabla1k_1e(double* f, double* g,
                            FINT li, FINT lj, FINT lk, CINTEnvVars envs)
        {
            FINT dj = envs.g_stride_j;
            FINT dk = envs.g_stride_k;
            double ak2 = -2 * envs.ak;
            FINT i, j, k, ptr;
            double* gx = g;
            double* gy = g + envs.g_size;
            double* gz = g + envs.g_size * 2;
            double* fx = f;
            double* fy = f + envs.g_size;
            double* fz = f + envs.g_size * 2;

            for (j = 0; j <= lj; j++)
            {
                ptr = dj * j;
                for (i = ptr; i <= ptr + li; i++)
                {
                    fx[i] = ak2 * gx[i + dk];
                    fy[i] = ak2 * gy[i + dk];
                    fz[i] = ak2 * gz[i + dk];
                }
            }
            for (k = 1; k <= lk; k++)
            {
                for (j = 0; j <= lj; j++)
                {
                    ptr = dj * j + dk * k;
                    for (i = ptr; i <= ptr + li; i++)
                    {
                        fx[i] = k * gx[i - dk] + ak2 * gx[i + dk];
                        fy[i] = k * gy[i - dk] + ak2 * gy[i + dk];
                        fz[i] = k * gz[i - dk] + ak2 * gz[i + dk];
                    }
                }
            }
        }


        /*
         * ( x^1 i j | k )
         * ri is the shift from the center R_O to the center of |i>
         * r - R_O = (r-R_i) + ri, ri = R_i - R_O
         */
        internal static void CINTx1i_1e(double* f, double* g, Vector3 ri,
                        FINT li, FINT lj, FINT lk, CINTEnvVars envs)
        {
            FINT i, j, k, ptr;
            FINT dj = envs.g_stride_j;
            FINT dk = envs.g_stride_k;
            double* gx = g;
            double* gy = g + envs.g_size;
            double* gz = g + envs.g_size * 2;
            double* fx = f;
            double* fy = f + envs.g_size;
            double* fz = f + envs.g_size * 2;

            for (k = 0; k <= lk; k++)
            {
                for (j = 0; j <= lj; j++)
                {
                    ptr = dj * j + dk * k;
                    for (i = ptr; i <= ptr + li; i++)
                    {
                        fx[i] = gx[i + 1] + ri[0] * gx[i];
                        fy[i] = gy[i + 1] + ri[1] * gy[i];
                        fz[i] = gz[i + 1] + ri[2] * gz[i];
                    }
                }
            }
        }

        internal static void CINTx1j_1e(double* f, double* g, Vector3 rj,
                        FINT li, FINT lj, FINT lk, CINTEnvVars envs)
        {
            FINT i, j, k, ptr;
            FINT dj = envs.g_stride_j;
            FINT dk = envs.g_stride_k;
            double* gx = g;
            double* gy = g + envs.g_size;
            double* gz = g + envs.g_size * 2;
            double* fx = f;
            double* fy = f + envs.g_size;
            double* fz = f + envs.g_size * 2;

            for (k = 0; k <= lk; k++)
            {
                for (j = 0; j <= lj; j++)
                {
                    ptr = dj * j + dk * k;
                    for (i = ptr; i <= ptr + li; i++)
                    {
                        fx[i] = gx[i + dj] + rj[0] * gx[i];
                        fy[i] = gy[i + dj] + rj[1] * gy[i];
                        fz[i] = gz[i + dj] + rj[2] * gz[i];
                    }
                }
            }
        }

        internal static void CINTx1k_1e(double* f, double* g, double* rk,
                        FINT li, FINT lj, FINT lk, CINTEnvVars envs)
        {
            FINT i, j, k, ptr;
            FINT dj = envs.g_stride_j;
            FINT dk = envs.g_stride_k;
            double* gx = g;
            double* gy = g + envs.g_size;
            double* gz = g + envs.g_size * 2;
            double* fx = f;
            double* fy = f + envs.g_size;
            double* fz = f + envs.g_size * 2;

            for (k = 0; k <= lk; k++)
            {
                for (j = 0; j <= lj; j++)
                {
                    ptr = dj * j + dk * k;
                    for (i = ptr; i <= ptr + li; i++)
                    {
                        fx[i] = gx[i + dk] + rk[0] * gx[i];
                        fy[i] = gy[i + dk] + rk[1] * gy[i];
                        fz[i] = gz[i + dk] + rk[2] * gz[i];
                    }
                }
            }
        }

        /*
         * gc    contracted GTO integral
         * nf    number of primitive integral
         * gp    primitive GTO integral
         * inc   increment of gp
         * shl   nth shell
         * ip    ith-1 primitive GTO
         */
        internal static void CINTprim_to_ctr(double* gc, FINT nf, double* gp,
                             FINT inc, FINT nprim, FINT nctr, double* coeff)
        {
            FINT n, i, k;
            double* pgc = gc;
            double c;

            for (i = 0; i < inc; i++)
            {
                //dger(nf, nctr, 1.d0, gp(i+1), inc, env(ptr), nprim, gc(1,i*nctr+1), nf)
                for (n = 0; n < nctr; n++)
                {
                    c = coeff[nprim * n];
                    if (c != 0)
                    {
                        for (k = 0; k < nf; k++)
                        {
                            pgc[k] += c * gp[k * inc + i];
                        }
                    }
                    // next cgto block
                    pgc += nf;
                }
            }
        }

        internal static void CINTprim_to_ctr_0(double* gc, double* gp, double* coeff, FINT nf,
                               FINT nprim, FINT nctr, FINT non0ctr, FINT* sortedidx)
        {
            FINT i;
            FINT n;
            double c0;

            for (i = 0; i < nctr; i++)
            {
                c0 = coeff[nprim * i];
                for (n = 0; n < nf; n++)
                {
                    gc[nf * i + n] = c0 * gp[n];
                }
            }
        }

        internal static void CINTprim_to_ctr_1(double* gc, double* gp, double* coeff, int nf,
                               FINT nprim, FINT nctr, FINT non0ctr, FINT* sortedidx)
        {
            FINT i, j;
            int n;
            double c0;

            for (i = 0; i < non0ctr; i++)
            {
                c0 = coeff[nprim * sortedidx[i]];
                j = sortedidx[i];
                for (n = 0; n < nf; n++)
                {
                    gc[nf * j + n] += c0 * gp[n];
                }
            }
        }

        /*
         * to optimize memory copy in cart2sph.c, remove the common factor for s
         * and p function in cart2sph
         */
        internal static double CINTcommon_fac_sp(FINT l)
        {
            switch (l)
            {
                case 0: return 0.282094791773878143;
                case 1: return 0.488602511902919921;
                default: return 1;
            }
        }


    }
}
