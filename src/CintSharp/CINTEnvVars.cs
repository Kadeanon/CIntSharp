using CintSharp.DataStructures;
using CintSharp.DataStructures.Native;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp
{
    internal unsafe class CINTEnvVars
    {
        internal Atm[] atm;
        internal Bas[] bas;
        internal double[] env;
        internal FINT[] shls;
        internal FINT natm;
        internal FINT nbas;

        internal FINT i_l;
        internal FINT j_l;
        internal FINT k_l;
        internal FINT l_l;
        internal FINT nfi;  // number of cartesian components
        internal FINT nfj;
        // in int1e_grids, the grids_offset and the number of grids
        internal FINT nfk;
        internal FINT grids_offset;
        internal FINT nfl;
        internal FINT ngrids;
        internal FINT nf;  // = nfi*nfj*nfk*nfl;
        internal FINT rys_order; // = nrys_roots for regular ERIs. can be nrys_roots/2 for SR ERIs
        internal FINT[] x_ctr;

        internal FINT gbits;
        internal FINT ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        internal FINT ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint.h
        internal FINT ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        internal FINT li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        internal FINT lj_ceil;
        internal FINT lk_ceil;
        internal FINT ll_ceil;
        internal FINT g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        internal FINT g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        internal FINT g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        internal FINT g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        internal FINT nrys_roots;
        internal FINT g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        internal FINT g2d_ijmax;
        internal FINT g2d_klmax;
        internal double common_factor;
        internal double expcutoff;
        internal Vector3 rirj; // diff by sign in different g0_2d4d algorithm
        internal Vector3 rkrl;
        internal Vector3 rx_in_rijrx;
        internal Vector3 rx_in_rklrx;

        internal Vector3 ri;
        internal Vector3 rj;
        internal Vector3 rk;
        // in int2e or int3c2e, the coordinates of the fourth shell
        // in int1e_grids, the pointer for the grids coordinates
        internal double[] rl;
        internal double[] grids;

        internal Func<FINT> f_g0_2e;
        internal Action f_g0_2d4d;
        internal Action<Span<double>, Span<double>, Span<FINT>, CINTEnvVars, bool> f_gout;
        internal CINTOpt[] opt;

        /* values are assigned during calculation */
        internal int[] idx;
        internal double ai;
        internal double aj;
        internal double ak;
        internal double al;
        internal double fac;
        internal Vector3 rij;
        internal Vector3 rkl;

        private CINTEnvVars()
        {
        }


        public static CINTEnvVars FromInt1e(Span<int> ng, int[] shls, Atm[] atm, int natm, Bas[] bas, int nbas, double[] env)
        {
            FINT i_sh = shls[0];
            FINT j_sh = shls[1];
            CINTEnvVars envs = new()
            {
                natm = natm,
                nbas = nbas,
                atm = atm,
                bas = bas,
                env = env,
                shls = shls,
                i_l = bas[i_sh].angleOf,
                j_l = bas[j_sh].angleOf,
                x_ctr = [bas[i_sh].numOfCont, bas[i_sh].numOfCont, 0, 0]
            };
            envs.nfi = (envs.i_l + 1) * (envs.i_l + 2) / 2;
            envs.nfj = (envs.j_l + 1) * (envs.j_l + 2) / 2;
            envs.nf = envs.nfi * envs.nfj;
            envs.common_factor = 1;
            if (env[PTR_EXPCUTOFF] == 0)
            {
                envs.expcutoff = EXPCUTOFF;
            }
            else
            {
                envs.expcutoff = Math.Max(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
            }

            envs.li_ceil = envs.i_l + ng[IINC];
            envs.lj_ceil = envs.j_l + ng[JINC];
            envs.ri = Unsafe.As<double, Vector3>(ref env[atm[bas[i_sh].atomOf].pointerForCoords]);
            envs.rj = Unsafe.As<double, Vector3>(ref env[atm[bas[j_sh].atomOf].pointerForCoords]);

            envs.gbits = ng[GSHIFT];
            envs.ncomp_e1 = ng[POS_E1];
            envs.ncomp_tensor = ng[TENSOR];
            if (ng[SLOT_RYS_ROOTS] > 0)
            {
                envs.nrys_roots = ng[SLOT_RYS_ROOTS];
            }
            else
            {
                envs.nrys_roots = (envs.li_ceil + envs.lj_ceil) / 2 + 1;
            }

            FINT dli, dlj;
            bool ibase = envs.li_ceil > envs.lj_ceil;
            if (ibase)
            {
                dli = envs.li_ceil + envs.lj_ceil + 1;
                dlj = envs.lj_ceil + 1;
                envs.rirj = envs.ri - envs.rj;
            }
            else
            {
                dli = envs.li_ceil + 1;
                dlj = envs.li_ceil + envs.lj_ceil + 1;
                envs.rirj = envs.rj - envs.ri;
            }
            envs.g_stride_i = envs.nrys_roots;
            envs.g_stride_j = envs.nrys_roots * dli;
            envs.g_size = envs.nrys_roots * dli * dlj;
            envs.g_stride_k = envs.g_size;
            envs.g_stride_l = envs.g_size;

            Debug.Assert(i_sh < SHLS_MAX);
            Debug.Assert(j_sh < SHLS_MAX);
            Debug.Assert(envs.i_l < ANG_MAX);
            Debug.Assert(envs.j_l < ANG_MAX);
            Debug.Assert(bas[i_sh].atomOf >= 0);
            Debug.Assert(bas[j_sh].atomOf >= 0);
            Debug.Assert(bas[i_sh].atomOf < natm);
            Debug.Assert(bas[j_sh].atomOf < natm);
            Debug.Assert(envs.nrys_roots < MXRYSROOTS);

            return envs;
        }

    }
}
