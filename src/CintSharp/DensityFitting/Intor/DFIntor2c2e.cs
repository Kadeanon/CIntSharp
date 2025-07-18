using CintSharp.Intor;
using CintSharp.Native.Libcint;
using SimpleHelpers.MultiAlg;
using System.Buffers;

namespace CintSharp.DensityFitting.Intor
{
    internal class DFIntor2c2e : DFIntorBase
    {
        public DFIntor2c2e(DFEnvs envs) : base(envs, "int2c2e")
        {
            if (type != IntorType.Spheric)
            {
                throw new ArgumentException($"The intor type of int2c2e is not spheric.");
            }
        }
        public override NDArray Invoke()
        {
            var intor = LibcintHandler.CreateIntor(Envs, $"{IntorName}_sph");
            var shellLengths = Envs.ShellLengths;
            var shellAuxLengths = Envs.ShellLengthsAux;
            var maxLength = shellAuxLengths.Max();
            var nshl = shellLengths.Length;
            var nshlAux = shellAuxLengths.Length;
            var ranges = Envs.RangesByShellsAux;
            var result = NDArray.CreateUninitialized([Components, Envs.NAOAux, Envs.NAOAux]);
            maxLength = maxLength * maxLength * Components;
            double[] caches = ArrayPool<double>.Shared.Rent(1024 * 512);
            double[] buffer = ArrayPool<double>.Shared.Rent(maxLength);
            int[] dims = ArrayPool<int>.Shared.Rent(3);
            int[] shls = ArrayPool<int>.Shared.Rent(3);
            for (int i = 0; i < nshlAux; i++)
            {
                int lengthI = shellAuxLengths[i];
                dims[0] = lengthI;
                shls[0] = i + nshl;
                for (int j = 0; j < nshlAux; j++)
                {
                    int lengthJ = shellAuxLengths[j];
                    dims[1] = lengthJ;
                    shls[1] = j + nshl;
                    intor.Invoke(buffer, dims, shls, Envs.Atms, Envs.Natm, Envs.Bases, Envs.Nbas, Envs.Envs, Optimizer, caches);
                    result[.., ranges[i], ranges[j]] =
                        new NDArray(buffer, [Components, lengthI, lengthJ],
                            [lengthI * lengthJ,
                             1,
                             lengthI]);
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
