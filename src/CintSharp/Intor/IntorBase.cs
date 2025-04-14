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
    /// <summary>
    /// A base class for all integrals.
    /// Inherit this class to create a kind of new integral.
    /// And implement the <see cref="Invoke"/> method to get the result.
    /// Check <seealso cref="Intor1e"/> and <seealso cref="Intor2e"/> for examples.
    /// </summary>
    public abstract class IntorBase : IDisposable
    {
        /// <summary>
        /// The <see cref="CIntEnvs"/> object that contains the params and some precompute properties.
        /// </summary>
        protected CIntEnvs Envs { get; }

        /// <summary>
        /// A flag to indicate whether the optimizer should be used.
        /// </summary>
        private bool ShouldOptimize { get; }

        /// <summary>
        /// The name of the integral without the suffix "_sph/_cart/spin". 
        /// It will be used to call the corresponding integral function in the libcint library.
        /// Add suffix to it manually to get the complete integrator name.
        /// </summary>
        protected string IntorName { get; }


        private nint optimizer;

        /// <summary>
        /// The optimizer pointer that can be used to optimize the integral calculation.
        /// Only be used when the <seealso cref="ShouldOptimize"/> is <see cref="true"/>.
        /// </summary>
        protected nint Optimizer => optimizer;

        /// <summary>
        /// The number of components of the integral.
        /// Used to determine the shape of the result tensor as gradient intor.
        /// </summary>
        protected int Components { get; }

        /// <summary>
        /// The type of the integral.
        /// For now, it has only been implemented for the spheric type.
        /// </summary>
        protected IntorType type;

        private bool disposedValue;

        /// <summary>
        /// Create a new integral object.
        /// </summary>
        /// <param name="envs">The <see cref="CIntEnvs"/> object to integrate.</param>
        /// <param name="intorName">The intor name, can</param>
        /// <param name="shouldOpt"></param>
        /// <exception cref="ArgumentException"></exception>
        protected IntorBase(CIntEnvs envs, string intorName, bool shouldOpt = true)
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

        /// <summary>
        /// Invoke the integral calculation.
        /// </summary>
        /// <returns>The <see cref="Tensor{double}"/> result that .</returns>
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
