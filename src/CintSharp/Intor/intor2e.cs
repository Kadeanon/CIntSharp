using CintSharp.DataStructures;
using CintSharp.Native.Libcint;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Intor
{
    public sealed class Intor2e : IntorBase
    {
        public Intor2e(CIntEnvs envs, string intorName) : base(envs, intorName)
        {
            if (type != IntorType.Spheric)
            {
                throw new ArgumentException($"The intor type of {intorName} is not spheric.");
            }
        }

        public override Tensor<double> Invoke()
        {
            var intor = LibcintHandler.CreateIntor(Envs, $"{IntorName}_sph");
            var shellLength = Envs.ShellLengths;
            int maxLength = shellLength.Max();
            var nshl = shellLength.Length;
            var ranges = Envs.RangesByShells;
            Tensor<double> result = Tensor.CreateUninitialized<double>([Components, Envs.NAO, Envs.NAO, Envs.NAO, Envs.NAO]);
            maxLength = maxLength * maxLength * maxLength * maxLength;
            maxLength *= Components;
            double[] caches = ArrayPool<double>.Shared.Rent(1024 * 512);
            double[] buffer = ArrayPool<double>.Shared.Rent(maxLength);
            int[] dims = ArrayPool<int>.Shared.Rent(4);
            int[] shls = ArrayPool<int>.Shared.Rent(4);
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
                    for (int k = 0; k < nshl; k++)
                    {
                        int lengthK = shellLength[k];
                        dims[2] = lengthK;
                        shls[2] = k;
                        for (int l = 0; l < nshl; l++)
                        {
                            int lengthL = shellLength[l];
                            dims[3] = lengthL;
                            shls[3] = l;
                            intor.Invoke(buffer, dims, shls, Envs.Atms, Envs.Natm, Envs.Bases, Envs.Nbas, Envs.Envs, Optimizer, caches);
                            result[.., ranges[i], ranges[j], ranges[k], ranges[l]] =
                                Tensor.Create(buffer, [Components, lengthL, lengthK, lengthJ, lengthI])
                                .PermuteDimensions([0, 4, 3, 2, 1]); // Ugly, maybe.
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
