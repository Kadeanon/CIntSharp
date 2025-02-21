using CintSharp.DataStructures.Native;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using static CintSharp.CintCommon;

namespace CintSharp.Intors.Int1e
{
    internal unsafe static class Int1e
    {
        /*
         * 1e GTO integral basic loop for < i|j>, no 1/r
         */
        internal static bool CINT1e_loop(scoped ref double gctr, CINTEnvVars envs, Span<double> cache, FINT int1e_type)
        {
            FINT[] shls = envs.shls;
            var bas = envs.bas;
            double[] env = envs.env;
            FINT i_sh = shls[0];
            FINT j_sh = shls[1];
            FINT i_ctr = envs.x_ctr[0];
            FINT j_ctr = envs.x_ctr[1];
            FINT i_prim = bas[i_sh].numOfPrim;
            FINT j_prim = bas[j_sh].numOfPrim;
            var ai = MemoryMarshal.CreateSpan(ref env[bas[i_sh].pointerOfExps], i_prim);
            var aj = MemoryMarshal.CreateSpan(ref env[bas[j_sh].pointerOfExps], j_prim);
            var ci = MemoryMarshal.CreateSpan(ref env[bas[i_sh].pointerOfCoeff], i_prim);
            var cj = MemoryMarshal.CreateSpan(ref env[bas[j_sh].pointerOfCoeff], i_prim);
            FINT n_comp = envs.ncomp_e1 * envs.ncomp_tensor;

            double expcutoff = envs.expcutoff;
            ref PairData pdata_ij = ref Unsafe.NullRef<PairData>();
            Misc.MALLOC_INSTACK(ref cache, out Span<double> log_maxci, i_prim + j_prim);//TODO : check this for alignment
            Misc.MALLOC_INSTACK(ref cache, out Span<PairData> pdata_base, i_prim * j_prim);
            var log_maxcj = log_maxci[i_prim..];
            log_maxci = log_maxci[..i_prim];
            CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
            CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
            if (CINTset_pairdata(pdata_base, ai, aj, envs.ri, envs.rj,
                         log_maxci, log_maxcj, envs.li_ceil, envs.lj_ceil,
                         i_prim, j_prim, envs.rirj * envs.rirj, expcutoff, env) != 0)
            {
                return false;
            }

            double fac1i, fac1j, expij;
            FINT ip, jp;
            bool empty = true;
            bool gempty = true;
            bool iempty = true;
            bool jempty = true;
            double* rij;
            Misc.MALLOC_INSTACK(ref cache, out Span<FINT> idx, envs.nf * 3);
            CINTg1e_index_xyz(idx, envs);

            Misc.MALLOC_INSTACK(ref cache, out Span<FINT> non0ctri, i_prim + j_prim + i_prim * i_ctr + j_prim * j_ctr);
            var non0ctrj = non0ctri.Slice(i_prim, j_prim);
            var non0idxi = non0ctrj.Slice(i_prim + j_prim, i_prim * i_ctr);
            var non0idxj = non0ctrj.Slice(i_prim + j_prim + i_prim * i_ctr, j_prim * j_ctr);
            non0ctri = non0ctri[..i_prim];
            CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
            CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);

            FINT nc = i_ctr * j_ctr;
            // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
            FINT leng = envs.g_size * 3 * ((1 << envs.gbits) + 1);
            FINT lenj = envs.nf * nc * n_comp; // gctrj
            FINT leni = envs.nf * i_ctr * n_comp; // gctri
            FINT len0 = envs.nf * n_comp; // gout
            FINT len = leng + lenj + leni + len0;
            ref double gout = ref Unsafe.NullRef<double>();
            ref double gctri = ref Unsafe.NullRef<double>();
            ref double gctrj = ref Unsafe.NullRef<double>();
            Misc.MALLOC_INSTACK(ref cache, out Span<double> g, len);  // must be allocated last in this function
            ref double g1 = ref MemoryMarshal.GetReference(cache);
            if (n_comp == 1)
            {
                gctrj = gctr;
            }
            else
            {
                gctrj = g1;
                g1 = ref Unsafe.Add(ref g1, lenj);
            }
            if (j_ctr == 1)
            {
                gctri = ref gctrj;
                iempty = jempty;
            }
            else
            {
                gctri = ref g1;
                g1 = ref Unsafe.Add(ref g1, leni);
            }
            if (i_ctr == 1)
            {
                gout = ref gctri;
                gempty = iempty;
            }
            else
            {
                gout = ref g1;
            }

            double common_factor = envs.common_factor
                    * CINTcommon_fac_sp(envs.i_l) * CINTcommon_fac_sp(envs.j_l);

            pdata_ij = ref pdata_base[0];
            for (jp = 0; jp < j_prim; jp++)
            {
                envs.aj = aj[jp];
                if (j_ctr == 1)
                {
                    fac1j = common_factor * cj[jp];
                }
                else
                {
                    fac1j = common_factor;
                    iempty = true;
                }
                for (ip = 0; ip < i_prim; ip++, pdata_ij = ref Unsafe.Add(ref pdata_ij, 1))
                {
                    if (pdata_ij.cceij > expcutoff)
                    {
                        continue;
                    }
                    envs.ai = ai[ip];
                    expij = pdata_ij.eij;
                    envs.rij = pdata_ij.rij;
                    if (i_ctr == 1)
                    {
                        fac1i = fac1j * ci[ip] * expij;
                    }
                    else
                    {
                        fac1i = fac1j * expij;
                    }
                    envs.fac = fac1i;
                    make_g1e_gout(ref gout, g, idx, envs, gempty, int1e_type);
                    //PRIM2CTR0(i, gout, envs.nf * n_comp);
                    if (i_ctr > 1)
                    {
                        if (iempty)
                        {
                            CINTprim_to_ctr_0(gctri, gout, ci[ip],
                                              envs.nf * n_comp, i_prim, i_ctr,
                                              non0ctri[ip],
                                              non0idxi[ip * i_ctr]);
                        }
                        else
                        {
                            CINTprim_to_ctr_1(gctri, gout, ci[ip],
                                              envs.nf * n_comp, i_prim, i_ctr,
                                              non0ctri[ip],
                                              non0idxi[ip * i_ctr]);
                        }
                    }
                    iempty = false;
                }
                if (!iempty)
                {
                    //PRIM2CTR0(j, gctri, envs.nf * i_ctr * n_comp);

                    if (j_ctr > 1)
                    {
                        if (jempty)
                        {
                            CINTprim_to_ctr_0(gctrj, gctri, cj[jp],
                                              envs.nf * i_ctr * n_comp, j_prim, j_ctr,
                                              non0ctri[jp],
                                              non0idxi[jp * j_ctr]);
                        }
                        else
                        {
                            CINTprim_to_ctr_1(gctrj, gctri, ci[jp],
                                              envs.nf * i_ctr * n_comp, j_prim, j_ctr,
                                              non0ctri[jp],
                                              non0idxi[jp * j_ctr]);
                        }
                    }
                    jempty = false;
                }
            }

            if (n_comp > 1 && !jempty)
            {
                CINTdmat_transpose(gctr, gctrj, envs.nf * nc, n_comp);
            }
            return !jempty;
        }


        public static CACHE_SIZE_T int1e_cache_size(CINTEnvVars envs)
        {
            FINT[] shls = envs.shls;
            var bas = envs.bas;
            FINT i_prim = bas[shls[0]].numOfPrim;
            FINT j_prim = bas[shls[1]].numOfPrim;
            FINT[] x_ctr = envs.x_ctr;
            FINT nc = envs.nf * x_ctr[0] * x_ctr[1];
            FINT n_comp = envs.ncomp_e1 * envs.ncomp_tensor;
            FINT leng = envs.g_size * 3 * ((1 << envs.gbits) + 1);
            FINT lenj = envs.nf * nc * n_comp;
            FINT leni = envs.nf * x_ctr[0] * n_comp;
            FINT len0 = envs.nf * n_comp;
            FINT pdata_size = i_prim * j_prim * 5
                               + i_prim * x_ctr[0]
                               + j_prim * x_ctr[1]
                               + (i_prim + j_prim) * 2 + envs.nf * 3;
            FINT cache_size = Math.Max(nc * n_comp + leng + lenj + leni + len0 + pdata_size,
                                 nc * n_comp + envs.nf * 8 * OF_CMPLX);
            return cache_size;
        }

        /*
         * 1e integrals <i|O|j> without 1/r
         */
        internal static int CINT1e_drv(Span<double> output, Span<FINT> dims, CINTEnvVars envs,
                       Span<double> cache, Action f_c2s, FINT int1e_type)
        {
            if (output.IsEmpty)
            {
                return int1e_cache_size(envs);
            }
            FINT[] x_ctr = envs.x_ctr;
            FINT nc = envs.nf * x_ctr[0] * x_ctr[1];
            FINT n_comp = envs.ncomp_e1 * envs.ncomp_tensor;
            double[]? stack = null;
            if (cache.IsEmpty)
            {
                int cache_size = int1e_cache_size(envs);
                stack = ArrayPool<double>.Shared.Rent(cache_size);
                cache = stack;
            }
            Misc.MALLOC_INSTACK(ref cache, out Span<double> gctr, nc * n_comp);

            bool has_value = CINT1e_loop(ref gctr[0], envs, cache, int1e_type);

            Span<FINT> counts = stackalloc FINT[4];
            if (dims.IsEmpty)
            {
                dims = counts;
            }
            if (f_c2s == c2s_sph_1e)
            {
                counts[0] = (envs.i_l * 2 + 1) * x_ctr[0];
                counts[1] = (envs.j_l * 2 + 1) * x_ctr[1];
            }
            else if (f_c2s == c2s_cart_1e)
            {
                counts[0] = envs.nfi * x_ctr[0];
                counts[1] = envs.nfj * x_ctr[1];
            }
            counts[2] = 1;
            counts[3] = 1;
            FINT nout = dims[0] * dims[1];
            FINT n;
            if (has_value)
            {
                for (n = 0; n < n_comp; n++)
                {
                    f_c2s(output[nout * n], gctr[nc * n], dims, envs, cache);
                }
            }
            else
            {
                for (n = 0; n < n_comp; n++)
                {
                    c2s_dset0(output[nout * n], dims, counts);
                }
            }
            if (stack != null)
            {
                ArrayPool<double>.Shared.Return(stack);
            }
            return has_value ? 1 : 0;
        }

        internal static CACHE_SIZE_T CINT1e_spinor_drv(Span<Complex> output, Span<FINT> dims, CINTEnvVars envs,
                               Span<double> cache, Action f_c2s, FINT int1e_type)
        {
            if (output.IsEmpty)
            {
                return int1e_cache_size(envs);
            }
            FINT[] x_ctr = envs.x_ctr;
            FINT nc = envs.nf * x_ctr[0] * x_ctr[1] * envs.ncomp_e1;
            double[]? stack = null;
            if (cache.IsEmpty)
            {
                int cache_size = int1e_cache_size(envs);
                stack = ArrayPool<double>.Shared.Rent(cache_size);
                cache = stack;
            }
            Misc.MALLOC_INSTACK(ref cache, out Span<double> gctr, nc * envs.ncomp_tensor);

            bool has_value = CINT1e_loop(ref gctr[0], envs, cache, int1e_type);

            Span<FINT> counts = stackalloc FINT[4];
            if (dims.IsEmpty)
            {
                dims = counts;
            }
            counts[0] = CINTcgto_spinor(envs.shls[0], envs.bas);
            counts[1] = CINTcgto_spinor(envs.shls[1], envs.bas);
            counts[2] = 1;
            counts[3] = 1;
            FINT nout = dims[0] * dims[1];
            FINT n;
            if (has_value)
            {
                for (n = 0; n < envs.ncomp_tensor; n++)
                {
                    f_c2s(output[nout * n], gctr[nc * n], dims, envs, cache);
                }
            }
            else
            {
                for (n = 0; n < envs.ncomp_tensor; n++)
                {
                    c2s_zset0(output[nout * n], dims, counts);
                }
            }

            if (stack != null)
            {
                ArrayPool<double>.Shared.Return(stack);
            }
            return has_value ? 1 : 0;
        }

        internal static void make_g1e_gout(Span<double> gout, Span<double> g, Span<FINT> idx,
                                  CINTEnvVars envs, bool empty, FINT int1e_type)
        {
            FINT ia;
            switch (int1e_type)
            {
                case 0:
                    G1e.CINTg1e_ovlp(g, envs);
                    envs.f_gout(gout, g, idx, envs, empty);
                    break;
                case 1:
                    G1e.CINTg1e_nuc(g, envs, -1);
                    envs.f_gout(gout, g, idx, envs, empty);
                    break;
                case 2:
                    for (ia = 0; ia < envs.natm; ia++)
                    {
                        G1e.CINTg1e_nuc(g, envs, ia);
                        envs.f_gout(gout, g, idx, envs, empty && ia == 0);
                    }
                    break;
            }
        }
        internal static void CINTgout1e(Span<double> gout, Span<double> g, Span<FINT> idx, CINTEnvVars envs, bool empty)
        {
            FINT nf = envs.nf;
            FINT n, ix, iy, iz;
            if (empty)
            {
                for (n = 0; n < nf; n++)
                {
                    ix = idx[n * 3 + 0];
                    iy = idx[n * 3 + 1];
                    iz = idx[n * 3 + 2];
                    gout[n] = g[ix] * g[iy] * g[iz];
                }
            }
            else
            {
                for (n = 0; n < nf; n++)
                {
                    ix = idx[n * 3 + 0];
                    iy = idx[n * 3 + 1];
                    iz = idx[n * 3 + 2];
                    gout[n] += g[ix] * g[iy] * g[iz];
                }
            }
        }

        internal static void CINTgout1e_nuc(Span<double> gout, Span<double> g, Span<FINT> idx, CINTEnvVars envs, bool empty)
        {
            FINT nf = envs.nf;
            FINT nrys_roots = envs.nrys_roots;
            FINT n, i;
            Span<double> gx, gy, gz;
            double s;

            if (empty)
            {
                for (n = 0; n < nf; n++)
                {
                    gx = g.Slice(idx[n * 3 + 0], nrys_roots);
                    gy = g.Slice(idx[n * 3 + 1], nrys_roots);
                    gz = g.Slice(idx[n * 3 + 2], nrys_roots);
                    s = 0;
                    for (i = 0; i < nrys_roots; i++)
                    {
                        s += gx[i] * gy[i] * gz[i];
                    }
                    gout[n] = s;
                }
            }
            else
            {
                for (n = 0; n < nf; n++)
                {
                    gx = g.Slice(idx[n * 3 + 0], nrys_roots);
                    gy = g.Slice(idx[n * 3 + 1], nrys_roots);
                    gz = g.Slice(idx[n * 3 + 2], nrys_roots);
                    s = 0;
                    for (i = 0; i < nrys_roots; i++)
                    {
                        s += gx[i] * gy[i] * gz[i];
                    }
                    gout[n] += s;
                }
            }
        }

        internal static CACHE_SIZE_T int1e_ovlp_sph(double[] output, FINT[] dims, FINT[] shls, Atm[] atm, FINT natm,
                             Bas[] bas, FINT nbas, double[] env, CINTOpt[] opt, double[] cache)
        {
            Span<int> ng = [0, 0, 0, 0, 0, 1, 1, 1];
            CINTEnvVars envs = CINTEnvVars.FromInt1e(ng, shls, atm, natm, bas, nbas, env);
            envs.f_gout = CINTgout1e;
            return CINT1e_drv(output, dims, envs, cache, c2s_sph_1e, 0);
        }

        internal static CACHE_SIZE_T int1e_ovlp_cart(double[] output, FINT[] dims, FINT[] shls, Atm[] atm, FINT natm,
                             Bas[] bas, FINT nbas, double[] env, CINTOpt[] opt, double[] cache)
        {
            Span<FINT> ng = [0, 0, 0, 0, 0, 1, 1, 1];
            CINTEnvVars envs = CINTEnvVars.FromInt1e(ng, shls, atm, natm, bas, nbas, env);
            envs.f_gout = CINTgout1e;
            return CINT1e_drv(output, dims, envs, cache, c2s_cart_1e, 0);
        }

        internal static CACHE_SIZE_T int1e_ovlp_spinor(Complex output, FINT[] dims, FINT[] shls, Atm[] atm, FINT natm,
                             Bas[] bas, FINT nbas, double[] env, CINTOpt[] opt, double[] cache)
        {
            Span<FINT> ng = [0, 0, 0, 0, 0, 1, 1, 1];
            CINTEnvVars envs = CINTEnvVars.FromInt1e(ng, shls, atm, natm, bas, nbas, env);
            envs.f_gout = CINTgout1e;
            return CINT1e_spinor_drv(output, dims, envs, cache, c2s_sf_1e, 0);
        }

        internal static void int1e_ovlp_optimizer(out CINTOpt? opt, FINT* atm, FINT natm,
                                  FINT* bas, FINT nbas, double* env)
        {
            opt = null;
        }

        internal static CACHE_SIZE_T int1e_nuc_sph(double[] output, FINT[] dims, FINT[] shls, Atm[] atm, FINT natm,
                             Bas[] bas, FINT nbas, double[] env, CINTOpt[] opt, double[] cache)
        {
            Span<FINT> ng = [0, 0, 0, 0, 0, 1, 0, 1];
            CINTEnvVars envs = CINTEnvVars.FromInt1e(ng, shls, atm, natm, bas, nbas, env);
            envs.f_gout = CINTgout1e_nuc;
            return CINT1e_drv(output, dims, envs, cache, c2s_sph_1e, 2);
        }

        internal static CACHE_SIZE_T int1e_nuc_cart(double[] output, FINT[] dims, FINT[] shls, Atm[] atm, FINT natm,
                             Bas[] bas, FINT nbas, double[] env, CINTOpt[] opt, double[] cache)
        {
            Span<FINT> ng = [0, 0, 0, 0, 0, 1, 0, 1];
            CINTEnvVars envs = CINTEnvVars.FromInt1e(ng, shls, atm, natm, bas, nbas, env);
            envs.f_gout = CINTgout1e_nuc;
            return CINT1e_drv(output, dims, &envs, cache, c2s_cart_1e, 2);
        }

        internal static CACHE_SIZE_T int1e_nuc_spinor(Complex[] output, FINT[] dims, FINT[] shls, Atm[] atm, FINT natm,
                             Bas[] bas, FINT nbas, double[] env, CINTOpt[] opt, double[] cache)
        {
            Span<FINT> ng = [0, 0, 0, 0, 0, 1, 0, 1];
            CINTEnvVars envs = CINTEnvVars.FromInt1e(ng, shls, atm, natm, bas, nbas, env);
            envs.f_gout = CINTgout1e_nuc;
            return CINT1e_spinor_drv(output, dims, envs, cache, c2s_sf_1e, 2);
        }

        internal static void int1e_nuc_optimizer(out CINTOpt? opt, FINT* atm, FINT natm,
                                  FINT* bas, FINT nbas, double* env)
        {
            opt = null;
        }

    }
}
