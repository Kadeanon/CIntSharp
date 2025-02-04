using CintSharp.DataStructures.Native;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DataStructures
{
    /// <summary>
    /// libcint中的shell
    /// </summary>
    public struct Shell
    {
        public int AngMomentum { get; }
        public int NumOfPrim { get; }
        public int NumOfCont { get; }
        public int ExpsPtr { get; }
        public int CoefsPtr { get; }

        public Shell(double[] exps, double[] coefs, int angMomentum, List<double> envs)
        {
            ExpsPtr = envs.Count;
            envs.AddRange(exps);
            CoefsPtr = envs.Count;
            envs.AddRange(coefs);
            NumOfPrim = exps.Length;
            NumOfCont = coefs.Length / exps.Length;
            AngMomentum = angMomentum;
        }

        public readonly void CopyToBasis(ref Bas bas, int atomIdx, List<double> envs)
        {
            bas.atomOf = atomIdx;
            bas.angleOf = AngMomentum;
            bas.numOfPrim = NumOfPrim;
            bas.numOfCont = NumOfCont;
            bas.pointerOfExps = ExpsPtr;
            bas.pointerOfCoeff = CoefsPtr;
        }

    }
}
