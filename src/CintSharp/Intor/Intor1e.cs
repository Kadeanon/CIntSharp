using CintSharp.DataStructures;
using CintSharp.Native.Libcint;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Intor
{
    public class Intor1e : IntorBase
    {
        public Intor1e(CIntEnvs envs, string intorName) : base(envs, intorName)
        {
        }

        public override Tensor<double> InvokeKernal(ReadOnlySpan<int> shellLength)
        {
            var intor = LibcintHandler.CreateIntor(Envs, $"{IntorName}_sph");
            int nshl = shellLength.Length;
            int[] shellOffset = Envs.OffsetsByShells;
            int maxLength = 0;
            Span<NRange> ranges = stackalloc NRange[nshl];
            for (int i = 0; i < nshl; i++)
            {
                ranges[i] = new NRange(shellOffset[i], shellOffset[i] + shellLength[i]);
                if (shellLength[i] > maxLength)
                {
                    maxLength = shellLength[i];
                }
            }
            Tensor<double> result = Tensor.CreateUninitialized<double>([Envs.NAO, Envs.NAO]);
            maxLength *= maxLength;
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
                    var resultChunk = result.AsTensorSpan(ranges[i], ranges[j]);
                    intor.Invoke(buffer, dims, shls, Envs.Atms, Envs.Natm, Envs.Bases, Envs.Nbas, Envs.Envs, Optimizer, caches);
                    for (int i2 = 0; i2 < lengthI; i2++)
                    {
                        for (int j2 = 0; j2 < lengthJ; j2++)
                        {
                            resultChunk[i2, j2] = buffer[i2 + j2 * lengthI];
                        }
                    }
                }
            }
            ArrayPool<double>.Shared.Return(caches);
            ArrayPool<double>.Shared.Return(buffer);
            ArrayPool<int>.Shared.Return(dims);
            ArrayPool<int>.Shared.Return(shls);
            return result;
        }
    }
}
