﻿using CintSharp.DataStructures;
using CintSharp.DataStructures.Native;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Native.Libcint
{
    public class LibcintHandler : DllHandler
    {
        public LibcintHandler(string path) : base(path)
        {
        }

        static LibcintHandler()
        {
            if (String.IsNullOrWhiteSpace(LibcintPath))
            {
                var attr = LibraryPathAttribute.Instance ?? throw new InvalidOperationException("The LibraryPathAttribute is not found in the assembly.");
                LibcintPath = attr.GetDllPlatformPath(Name);
            }
            Instance ??= new LibcintHandler(LibcintPath);
        }

        public static string Name => "Libcint";

        public static LibcintHandler Instance { get; private set; }

        delegate double CintGtoNormaler(int angMomentum, double v);

        public delegate void OptimizerCreator(ref IntPtr opt, Atm[] atm, int natm, Bas[] bas, int nbas, double[] env);

        public delegate void OptimizerDestroyer(ref IntPtr opt);

        public delegate int Intor(double[] output, int[] dims, int[] shls, Atm[] atm, int natm, Bas[] bas, int nbas, double[] env, IntPtr opt, double[] cache);

        CintGtoNormaler? CINTgto_normDelegate;

        public static double CINTgto_norm(int angMomentum, double exp)
        {
            if(Instance == null)
            {
                throw new InvalidOperationException("The LibcintHandler is not initialized.");
            }

            Instance.CINTgto_normDelegate ??= Instance.Invoke<CintGtoNormaler>("CINTgto_norm");

            return Instance.CINTgto_normDelegate(angMomentum, exp);
        }

        public static void GetOptimizer(ref IntPtr opt, CIntEnvs envs, string apiName)
        {
            if (Instance == null)
            {
                throw new InvalidOperationException("The LibcintHandler is not initialized.");
            }
            try
            {
                Instance.Invoke<OptimizerCreator>(apiName)(ref opt, envs.Atms, envs.Natm, envs.Bases, envs.Nbas, envs.Envs);
                return;
            }
            catch (Exception e) 
            {
                throw new InvalidOperationException($"Failed to create optimizer for {apiName}: {e.Message}", e);
            }
        }

        public static void ReleaseOptimizer(ref IntPtr opt) 
        {

            if (Instance == null)
            {
                throw new InvalidOperationException("The LibcintHandler is not initialized.");
            }
            try
            {
                Instance.Invoke<OptimizerDestroyer>("CINTdel_optimizer")(ref opt);
            }
            catch (Exception e)
            {
                throw new InvalidOperationException($"Failed to destroy optimizer: {e.Message}", e);
            }
        }

        public static Intor CreateIntor(CIntEnvs envs, string apiName)
        {
            if (Instance == null)
            {
                throw new InvalidOperationException("The LibcintHandler is not initialized.");
            }
            try
            {
                return Instance.Invoke<Intor>(apiName); 
            }
            catch (Exception e)
            {
                throw new InvalidOperationException($"Failed to create intor for {apiName}: {e.Message}", e);
            }
        }
    }
}
