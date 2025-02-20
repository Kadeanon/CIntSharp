global using FINT = int;
global using CACHE_SIZE_T = int;
global using static CintSharp.CintCommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CintSharp.DataStructures;
using CintSharp.DataStructures.Native;

namespace CintSharp
{
    public static class CintCommon
    {
        #region const
        public static readonly Version CINT_VERSION = new(6, 1, 1);
        public const int CINT_SOVERSION = 6;

        // global parameters in env
        // Overall cutoff for integral prescreening, value needs to be ~ln(threshold)
        internal const int PTR_EXPCUTOFF = 0;
        // R_C of (r-R_C) in dipole, GIAO operators
        internal const int PTR_COMMON_ORIG = 1;
        // R_O in 1/|r-R_O|
        internal const int PTR_RINV_ORIG = 4;
        // ZETA parameter for Gaussian charge distribution (Gaussian nuclear model)
        internal const int PTR_RINV_ZETA = 7;
        // omega parameter in range-separated coulomb operator
        // LR interaction: erf(omega*r12)/r12 if omega > 0
        // SR interaction: erfc(omega*r12)/r12 if omega < 0
        internal const int PTR_RANGE_OMEGA = 8;
        // Yukawa potential and Slater-type geminal e^{-zeta r}
        internal const int PTR_F12_ZETA = 9;
        // Gaussian type geminal e^{-zeta r^2}
        internal const int PTR_GTG_ZETA = 10;
        internal const int NGRIDS = 11;
        internal const int PTR_GRIDS = 12;
        internal const int PTR_ENV_START = 20;


        // slots of atm
        internal const int CHARGE_OF = 0;
        internal const int PTR_COORD = 1;
        internal const int NUC_MOD_OF = 2;
        internal const int PTR_ZETA = 3;
        internal const int PTR_FRAC_CHARGE = 4;
        internal const int RESERVE_ATMSLOT = 5;
        internal const int ATM_SLOTS = 6;


        // slots of bas
        internal const int ATOM_OF = 0;
        internal const int ANG_OF = 1;
        internal const int NPRIM_OF = 2;
        internal const int NCTR_OF = 3;
        internal const int KAPPA_OF = 4;
        internal const int PTR_EXP = 5;
        internal const int PTR_COEFF = 6;
        internal const int RESERVE_BASLOT = 7;
        internal const int BAS_SLOTS = 8;

        // slots of gout
        internal const int POSX = 0;
        internal const int POSY = 1;
        internal const int POSZ = 2;
        internal const int POS1 = 3;
        // For 2-electron integral with two spin operators
        // SIGMA1X * SIGMA2X     0
        // SIGMA1Y * SIGMA2X     1
        // SIGMA1Z * SIGMA2X     2
        // I1_2x2  * SIGMA2X     3
        // SIGMA1X * SIGMA2Y     4
        // SIGMA1Y * SIGMA2Y     5
        // SIGMA1Z * SIGMA2Y     6
        // I1_2x2  * SIGMA2Y     7
        // SIGMA1X * SIGMA2Z     8
        // SIGMA1Y * SIGMA2Z     9
        // SIGMA1Z * SIGMA2Z     10
        // I1_2x2  * SIGMA2Z     11
        // SIGMA1X * I2_2x2      12
        // SIGMA1Y * I2_2x2      13
        // SIGMA1Z * I2_2x2      14
        // I1_2x2  * I2_2x2      15
        internal const int POSXX = 0;
        internal const int POSYX = 1;
        internal const int POSZX = 2;
        internal const int POS1X = 3;
        internal const int POSXY = 4;
        internal const int POSYY = 5;
        internal const int POSZY = 6;
        internal const int POS1Y = 7;
        internal const int POSXZ = 8;
        internal const int POSYZ = 9;
        internal const int POSZZ = 10;
        internal const int POS1Z = 11;
        internal const int POSX1 = 12;
        internal const int POSY1 = 13;
        internal const int POSZ1 = 14;
        internal const int POS11 = 15;

        // tensor
        internal const int TSRX = 0;
        internal const int TSRY = 1;
        internal const int TSRZ = 2;
        internal const int TSRXX = 0;
        internal const int TSRXY = 1;
        internal const int TSRXZ = 2;
        internal const int TSRYX = 3;
        internal const int TSRYY = 4;
        internal const int TSRYZ = 5;
        internal const int TSRZX = 6;
        internal const int TSRZY = 7;
        internal const int TSRZZ = 8;

        // other boundaries
        internal const int MXRYSROOTS = 32; // > ANG_MAX*2+1 for 4c2e
        internal const int ANG_MAX = 15; // l = 0..15
        internal const int LMAX1 = 16; // > ANG_MAX
        internal const int CART_MAX = 136; // > (ANG_MAX*(ANG_MAX+1)/2)
        internal const int SHLS_MAX = 1048576;
        internal const int NPRIM_MAX = 64;
        internal const int NCTR_MAX = 64;

        internal const int POINT_NUC = 1;
        internal const int GAUSSIAN_NUC = 2;
        internal const int FRAC_CHARGE_NUC = 3;

        #endregion

        #region config
        internal const int IINC = 0;
        internal const int JINC = 1;
        internal const int KINC = 2;
        internal const int LINC = 3;
        internal const int GSHIFT = 4;
        internal const int POS_E1 = 5;
        internal const int POS_E2 = 6;
        internal const int SLOT_RYS_ROOTS = 6;
        internal const int TENSOR = 7;
        internal const int EXPCUTOFF = 60;
        // ~ 1e-15
        internal const int MIN_EXPCUTOFF = 40;
        internal const int OF_CMPLX = 2;
        internal const int GRID_BLKSIZE = 104;
        #endregion

        internal struct PairData
        {
            public Vector3 rij;
            public double eij;
            public double cceij;
        }
        internal unsafe class CINTOpt
        {
            FINT** index_xyz_array; // LMAX1**4 pointers to index_xyz
            FINT** non0ctr;
            FINT** sortedidx;
            FINT nbas;
            double** log_max_coeff;
            PairData** pairdata;  // NULL indicates not-initialized, NO_VALUE can be skipped
        }

        internal static void CINTcart_comp(Span<FINT> nx, Span<FINT> ny, Span<FINT> nz, FINT lmax)
{
        FINT inc = 0;
        FINT lx, ly, lz;

        for (lx = lmax; lx >= 0; lx--) {
                for (ly = lmax - lx; ly >= 0; ly--) {
                        lz = lmax - lx - ly;
                        nx[inc] = lx;
                        ny[inc] = ly;
                        nz[inc] = lz;
                        inc++;
                }
}
}
    }
}
