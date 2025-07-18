using CintSharp.DataStructures;
using SimpleHelpers;
using SimpleHelpers.LinearAlg;
using SimpleHelpers.MultiAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace CintSharp.DensityFitting.Intor
{

    public static partial class DFIntorUtils
    {
        #region Intors
        public static NDArray InvokeIntor(this DFEnvs envs, string intorName)
        {
            using var intor = envs.CreateIntor(intorName);
            return intor.Invoke();
        }

        public static NDArray GetOvlp(this DFEnvs envs) => envs.InvokeIntor("int1e_ovlp");

        public static NDArray GetKin(this DFEnvs envs) => envs.InvokeIntor("int1e_kin");

        public static NDArray GetNuc(this DFEnvs envs) => envs.InvokeIntor("int1e_nuc");

        public static NDArray GetERI(this DFEnvs envs) => envs.InvokeIntor("int2e");

        public static Matrix Get2c2e(this DFEnvs envs) => envs.InvokeIntor("int2c2e").AsMatrix();

        public static NDArray Get3c2e(this DFEnvs envs) => envs.InvokeIntor("int3c2e");

        public static NDArray GetHCore(this DFEnvs envs) =>
            envs.GetKin() + envs.GetNuc();

        #endregion
    }
}
