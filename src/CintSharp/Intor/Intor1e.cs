using CintSharp.DataStructures;
using CintSharp.Native.Libcint;
using System;
using System.Buffers;
using System.Collections.Generic;
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
            if (type != IntorUtils.IntorType.Spheric)
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
            Tensor<double> result = Tensor.CreateUninitialized<double>([Envs.NAO, Envs.NAO, Components]);
            maxLength *= maxLength;
            maxLength *= Components;
            double[] caches = ArrayPool<double>.Shared.Rent(1024);
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
                    var resultChunk = result.AsTensorSpan(ranges[i], ranges[j], ..);
                    intor.Invoke(buffer, dims, shls, Envs.Atms, Envs.Natm, Envs.Bases, Envs.Nbas, Envs.Envs, Optimizer, caches);
                    for (int i2 = 0; i2 < lengthI; i2++)
                    {
                        for (int j2 = 0; j2 < lengthJ; j2++)
                        {
                            for(int c = 0; c < Components; c++)
                            {
                                resultChunk[i2, j2, c] = buffer[i2 + j2 * lengthI];
                            }
                        }
                    }
                }
            }
            ArrayPool<int>.Shared.Return(shls);
            ArrayPool<int>.Shared.Return(dims);
            ArrayPool<double>.Shared.Return(caches);
            ArrayPool<double>.Shared.Return(buffer);
            if(Components == 1) 
            {
                result = result.Reshape([Envs.NAO, Envs.NAO]);
            }
            return result;
        }
    }
}
