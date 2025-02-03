using CintSharp.DataStructures;
using CIntSharp.DataStructures.Native;
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

        public delegate double CINTgto_normaler(int angMomentum, double v);

        public delegate void CintOptimizer(out IntPtr opt, Atm[] atm, int natm, Bas[] bas, int nbas, double[] env);

        public delegate int CintIntor(double[] output, int[] dims, int[] shls, Atm[] atm, int natm, Bas[] bas, int nbas, double[] env, IntPtr opt, double[] cache);

        CINTgto_normaler? CINTgto_normDelegate;

        public static double CINTgto_norm(int angMomentum, double exp)
        {
            if(Instance == null)
            {
                throw new InvalidOperationException("The LibcintHandler is not initialized.");
            }

            Instance.CINTgto_normDelegate ??= Instance.Invoke<CINTgto_normaler>("CINTgto_norm");

            return Instance.CINTgto_normDelegate(angMomentum, exp);
        }

        public static IntPtr GetOptimizer(CIntEnvs envs, string apiName)
        {
            if (Instance == null)
            {
                throw new InvalidOperationException("The LibcintHandler is not initialized.");
            }
            try
            {
                Instance.Invoke<CintOptimizer>(apiName)(out nint opt, envs.Atms, envs.Natm, envs.Bases, envs.Nbas, envs.Envs);
                return opt;
            }
            catch (Exception e) 
            {
                throw new InvalidOperationException($"Failed to create optimizer for {apiName}: {e.Message}", e);
            }
        }


        public static CintIntor CreateIntor(CIntEnvs envs, string apiName)
        {
            if (Instance == null)
            {
                throw new InvalidOperationException("The LibcintHandler is not initialized.");
            }
            CintIntor calculator;
            try
            {
                return Instance.Invoke<CintIntor>(apiName); 
            }
            catch (Exception e)
            {
                throw new InvalidOperationException($"Failed to create optimizer for {apiName}: {e.Message}", e);
            }
        }
    }
}
