using CintSharp.Intor;
using CintSharp.Native.Libcint;
using SimpleHelpers.MultiAlg;
using System.Buffers;

namespace CintSharp.DensityFitting.Intor
{
    internal class DFIntor3c2e : DFIntorBase
    {
        public DFIntor3c2e(DFEnvs envs) : base(envs, "int3c2e")
        {
            if (type != IntorType.Spheric)
            {
                throw new ArgumentException($"The intor type of int3c2e is not spheric.");
            }
        }
        public override NDArray Invoke()
        {
            var intor = LibcintHandler.CreateIntor(Envs, $"{IntorName}_sph");
            var shellLengths = Envs.ShellLengths;
            var shellAuxLengths = Envs.ShellLengthsAux;
            var maxLength = shellLengths.Max();
            var maxLengthAux = shellAuxLengths.Max();
            var nshl = shellLengths.Length;
            var nshlAux = shellAuxLengths.Length;
            var ranges = Envs.RangesByShells;
            var rangesAux = Envs.RangesByShellsAux;
            var result = NDArray.CreateUninitialized([Components, Envs.NAO, Envs.NAO, Envs.NAOAux]);
            maxLength *= maxLength * maxLengthAux * Components;
            double[] caches = ArrayPool<double>.Shared.Rent(1024 * 512);
            double[] buffer = ArrayPool<double>.Shared.Rent(maxLength);
            int[] dims = ArrayPool<int>.Shared.Rent(3);
            int[] shls = ArrayPool<int>.Shared.Rent(3);
            for (int i = 0; i < nshl; i++)
            {
                int lengthI = shellLengths[i];
                dims[0] = lengthI;
                shls[0] = i;
                for (int j = 0; j < nshl; j++)
                {
                    int lengthJ = shellLengths[j];
                    dims[1] = lengthJ;
                    shls[1] = j;
                    for (int k = 0; k < nshlAux; k++)
                    {
                        int lengthK = shellAuxLengths[k];
                        dims[2] = lengthK;
                        shls[2] = k + nshl;
                            intor.Invoke(buffer, dims, shls, Envs.Atms, Envs.Natm, Envs.Bases, Envs.Nbas, Envs.Envs, Optimizer, caches);
                            result[.., ranges[i], ranges[j], rangesAux[k]] =
                                new NDArray(buffer, [Components, lengthI, lengthJ, lengthK],
                                    [lengthI * lengthJ * lengthK,
                                     1,
                                     lengthI,
                                     lengthJ * lengthI]);
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
