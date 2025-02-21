using CintSharp.DataStructures;
using CintSharp.Native.Libcint;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Intor
{
    public abstract class IntorBase : IDisposable
    {
        protected CIntEnvs Envs { get; }
        private bool ShouldOptimize { get; }

        protected string IntorName { get; }
        private nint optimizer;

        protected nint Optimizer => optimizer;

        protected int Components { get; }

        protected IntorUtils.IntorType type;

        public IntorBase(CIntEnvs envs, string intorName, bool shouldOpt = true)
        {
            Components = IntorUtils.GetIntorComp(ref intorName, out type);
            if(Components == 0)
            {
                throw new ArgumentException($"The intor name {intorName} is not supported.");
            }
            ShouldOptimize = shouldOpt;
            Envs = envs;
            IntorName = intorName;
            if (ShouldOptimize)
            {
                LibcintHandler.GetOptimizer(ref optimizer, Envs, $"{IntorName}_optimizer");
            }
        }

        public Tensor<double> Invoke()
        {
            int nao = Envs.NAO;
            int nshl = Envs.OffsetsByShells.Length - 1;
            int[] shellLength = ArrayPool<int>.Shared.Rent(nshl);
            int[] shellOffset = Envs.OffsetsByShells;
            int lastOffset = Envs.OffsetsByShells[0];
            int currentOffset;
            for (int i = 0; i < nshl; i++)
            {
                currentOffset = Envs.OffsetsByShells[i + 1];
                int length = currentOffset - lastOffset;
                lastOffset = currentOffset;
                shellLength[i] = length;
            }
            Tensor<double> result = InvokeKernal(shellLength.AsSpan(0, nshl));
            ArrayPool<int>.Shared.Return(shellLength);
            return result;
        }

        protected abstract Tensor<double> InvokeKernal(ReadOnlySpan<int> shellLength);

        public void Dispose()
        {
            if(ShouldOptimize && optimizer != nint.Zero) 
            {
                LibcintHandler.ReleaseOptimizer(ref optimizer);
                optimizer = nint.Zero;
            }
            GC.SuppressFinalize(this);
        }
    }
}
