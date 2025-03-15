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
        private bool disposedValue;

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

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                //if (disposing)
                //{
                //    // TODO: 释放托管状态(托管对象)
                //}

                // TODO: 释放未托管的资源(未托管的对象)并重写终结器
                // TODO: 将大型字段设置为 null
                if (ShouldOptimize && optimizer != nint.Zero)
                {
                    LibcintHandler.ReleaseOptimizer(ref optimizer);
                }
                disposedValue = true;
            }
        }

        // TODO: 仅当“Dispose(bool disposing)”拥有用于释放未托管资源的代码时才替代终结器
        ~IntorBase()
        {
            // 不要更改此代码。请将清理代码放入“Dispose(bool disposing)”方法中
            Dispose(disposing: false);
        }

        public void Dispose()
        {
            // 不要更改此代码。请将清理代码放入“Dispose(bool disposing)”方法中
            Dispose(disposing: true);
            GC.SuppressFinalize(this);
        }
    }
}
