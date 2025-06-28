using CintSharp.DataStructures;
using CintSharp.Native.Libcint;
using SimpleHelpers.MultiAlg;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Intor
{
    public sealed class Intor1e : IntorBase
    {
        public Intor1e(CIntEnvs envs, string intorName) : base(envs, intorName)
        {
            if (type != IntorType.Spheric)
            {
                throw new ArgumentException($"The intor type of {intorName} is not spheric.");
            }
        }

        public override NDArray Invoke()
        {
            var intor = LibcintHandler.CreateIntor(Envs, $"{IntorName}_sph");
            var shellLength = Envs.ShellLengths;
            int maxLength = shellLength.Max();
            var nshl = shellLength.Length;
            var ranges = Envs.RangesByShells;
            NDArray result = NDArray.CreateUninitialized
                ([Components, Envs.NAO, Envs.NAO]);
            Tensor<double> val = Tensor.CreateUninitialized<double>
                ([Components, Envs.NAO, Envs.NAO]);
            maxLength *= maxLength;
            maxLength *= Components;
            double[] caches = ArrayPool<double>.Shared.Rent(1024 * Components * Components);
            double[] buffer = ArrayPool<double>.Shared.Rent(maxLength);
            int[] dims = ArrayPool<int>.Shared.Rent(2);
            int[] shls = ArrayPool<int>.Shared.Rent(2);
            for (int i = 0; i < nshl; i++)
            {
                int lengthI = shellLength[i];
                dims[0] = lengthI;
                shls[0] = i;
                for (int j = 0; j < nshl; j++)
                {
                    int lengthJ = shellLength[j];
                    dims[1] = lengthJ;
                    shls[1] = j;
                    intor.Invoke(buffer, dims, shls, 
                        Envs.Atms, Envs.Natm, Envs.Bases, Envs.Nbas, Envs.Envs, Optimizer, caches);

                    var subVal = val[.., ranges[i], ranges[j]];
                    var subResult = result[.., ranges[i], ranges[j]];

                    for (int c = 0; c < Components; c++)
                    {
                        for (int a = 0; a < lengthI; a++)
                        {
                            for (int b = 0; b < lengthJ; b++)
                            {
                                var valval = subVal[c, a, b];
                                var res = subResult[c, a, b];
                                Debug.Assert(Math.Abs(valval - res) < 1e-10,
                                    $"The value {val} and result {res} are not equal.");
                            }
                        }
                    }

                    val[.., ranges[i], ranges[j]] =
                        Tensor.Create(buffer, [Components, lengthJ, lengthI])
                        .PermuteDimensions([0, 2, 1]);
                    result[.., ranges[i], ranges[j]] =
                        new NDArray(buffer, [Components, lengthI, lengthJ],
                        [lengthI * lengthJ,
                        1,
                        lengthI]);

                    for (int c = 0; c < Components; c++)
                    {
                        for (int a = 0; a < lengthI; a++)
                        {
                            for (int b = 0; b < lengthJ; b++)
                            {
                                var valval = subVal[c, a, b];
                                var res = subResult[c, a, b];
                                Debug.Assert(Math.Abs(valval - res) < 1e-10,
                                    $"The value {val} and result {res} are not equal.");
                            }
                        }
                    }
                }
            }

            ArrayPool<int>.Shared.Return(shls);
            ArrayPool<int>.Shared.Return(dims);
            ArrayPool<double>.Shared.Return(caches);
            ArrayPool<double>.Shared.Return(buffer);
            if (Components == 1) 
            {
                result = result.SqueezeDimension(0);
            }
            return result;
        }
    }
}
