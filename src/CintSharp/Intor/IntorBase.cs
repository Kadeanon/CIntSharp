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

        public abstract Tensor<double> Invoke();

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
