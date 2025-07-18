using CintSharp.BasisParser;
using CintSharp.DataStructures;
using CintSharp.DensityFitting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Examples.LowScale.DensityFitting
{
    public static class DFUtils
    {
        public static DFBuilder WithDF(this EnvBuilder env,
            string basisName, IBasisParser? basisParser = null)
        {
            var dfBuilder = new DFBuilder(env, basisParser);
            return dfBuilder.SetBasis(basisName);
        }

        public static DFBuilder WithDF(this EnvBuilder env,
            Func<Atom, string> basisNameSetter, IBasisParser? basisParser = null)
        {
            var dfBuilder = new DFBuilder(env, basisParser);
            return dfBuilder.SetBasis(basisNameSetter);
        }
    }
}
