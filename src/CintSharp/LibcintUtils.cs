using CintSharp.DataStructures.Native;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp
{
    public static class LibcintUtils
    {
        public static int LenCart(int l)
        {
            return (l + 1) * (l + 2) / 2;
        }

        public static int LenSpinor(this Bas bas)
        {
            if (bas.kappaForSpinor == 0)
            {
                return 4 * bas.angleOf + 2;
            }
            else if (bas.kappaForSpinor < 0)
            {
                return 2 * bas.angleOf + 2;
            }
            else
            {
                return 2 * bas.angleOf;
            }
        }

        /* 
         * Num. of contracted cartesian GTO = 2j+1 * n_contraction
         */
        public static int CgtoCart(this Bas bas)
        {
            int l = bas.angleOf;
            return (l + 1) * (l + 2) / 2 * bas.numOfCont;
        }

        /* 
         * Num. of contracted spheric GTO = 2j+1 * n_contraction
         */
        public static int CgtoSpheric(this Bas bas)
        {
            return (bas.angleOf * 2 + 1) * bas.numOfCont;
        }

        /* 
         * Num. of contracted spinor GTO
         */
        public static int CgtoSpinor(this Bas bas)
        {
            return LenSpinor(bas) * bas.numOfCont;
        }

        /*
         * tot. primitive atomic spheric GTOs in a shell
         */
        public static int TotPgtoSpheric(Bas[] bases)
        {
            int s = 0;
            foreach(var bas in bases)
            {
                s += (bas.angleOf * 2 + 1) * bas.numOfPrim;
            }
            return s;
        }

        /*
         * tot. primitive atomic spinors in a shell
         */
        public static int TotPgtoSpinor(Bas[] bases)
        {
            int s = 0;

            foreach (var bas in bases)
            {
                s += LenSpinor(bas) * bas.numOfPrim;
            }
            return s;
        }

        public static int TotCgtoAccum(Func<Bas, int> f, Bas[] bases)
        {
            int s = 0;
            foreach (var bas in bases)
            {
                s += f(bas);
            }
            return s;
        }
        /*
         * tot. contracted atomic spheric GTOs in a shell
         */
        public static int TotCgtoSpheric(Bas[] bases)
        {
            return TotCgtoAccum(CgtoSpheric, bases);
        }

        /*
         * tot. contracted atomic spinors in a shell
         */
        public static int TotCgtoSpinor(Bas[] bases)
        {
            return TotCgtoAccum(CgtoSpinor, bases);
        }

        /*
         * tot. contracted atomic spinors in a shell
         */
        public static int TotCgtoCart(Bas[] bases)
        {
            return TotCgtoAccum(CgtoCart, bases);
        }

        public static void ShellsCgtoOffset(Func<Bas, int> f, int[] aoLoc,
                                       Bas[] bases)
        {
            int i;
            int nbas = bases.Length;
            aoLoc[0] = 0;
            for (i = 1; i < nbas; i++)
            {
                aoLoc[i] = aoLoc[i - 1] + f(bases[i - 1]);
            }
        }
        /*
         * offset of each shell for real spheric GTOs
         */
        public static void ShellsCartOffset(int[] aoLoc, Bas[] bases)
        {
            ShellsCgtoOffset(CgtoCart, aoLoc, bases);
        }

        /*
         * offset of each shell for real spheric GTOs
         */
        public static void ShellsSphericOffset(int[] aoLoc, Bas[] bases)
        {
            ShellsCgtoOffset(CgtoSpheric, aoLoc, bases);
        }

        /*
         * offset of each shell for AO spinors
         */
        public static void ShellsSpinorOffset(int[] aoLoc, Bas[] bases)
        {
            ShellsCgtoOffset(CgtoSpinor, aoLoc, bases);
        }


        /*
         * GTO = x^{nx}y^{ny}z^{nz}e^{-ar^2}
         */
        public static void CartComp(int[] nx, int[] ny, int[] nz, int lmax)
        {
            int inc = 0;
            int lx, ly, lz;

            for (lx = lmax; lx >= 0; lx--)
            {
                for (ly = lmax - lx; ly >= 0; ly--)
                {
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
