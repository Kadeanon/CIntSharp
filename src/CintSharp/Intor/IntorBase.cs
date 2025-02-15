﻿using CintSharp.DataStructures;
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
    public abstract class IntorBase : IDisposable
    {
        protected CIntEnvs Envs { get; }

        protected string IntorName { get; }

        protected nint Optimizer { get; private set; }

        protected int Components { get; }

        protected IntorUtils.IntorType type;

        public IntorBase(CIntEnvs envs, string intorName)
        {
            Components = IntorUtils.GetIntorComp(ref intorName, out type);
            if(Components == 0)
            {
                throw new ArgumentException($"The intor name {intorName} is not supported.");
            }
            Envs = envs;
            IntorName = intorName;
            Optimizer = LibcintHandler.GetOptimizer(Envs, $"{IntorName}_optimizer");
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
            var result = InvokeKernal(shellLength.AsSpan(0, nshl));
            ArrayPool<int>.Shared.Return(shellLength);
            return result;
        }

        protected abstract Tensor<double> InvokeKernal(ReadOnlySpan<int> shellLength);

        public void Dispose()
        {
            if(Optimizer != nint.Zero) 
            {
                LibcintHandler.ReleaseOptimizer(Optimizer);
                Optimizer = nint.Zero;
            }
        }
    }
}
