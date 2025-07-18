using CintSharp.DataStructures;
using CintSharp.Intor;
using CintSharp.Native.Libcint;
using SimpleHelpers.MultiAlg;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DensityFitting.Intor
{
    /// <summary>
    /// A base class for all integrals.
    /// Inherit this class to create a kind of new integral.
    /// And implement the <see cref="Invoke"/> method to get the result.
    /// Check <seealso cref="Intor1e"/> and <seealso cref="Intor2e"/> for examples.
    /// </summary>
    public abstract class DFIntorBase : IDisposable
    {
        /// <summary>
        /// The <see cref="DFEnvs"/> object that contains the params and some precompute properties.
        /// </summary>
        protected DFEnvs Envs { get; }

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
        protected DFIntorBase(DFEnvs envs, string intorName, bool shouldOpt = true)
        {
            Components = IntorUtils.GetIntorComp(ref intorName, out type);
            if (Components == 0)
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
        public abstract NDArray Invoke();

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (ShouldOptimize && optimizer != nint.Zero)
                {
                    LibcintHandler.ReleaseOptimizer(ref optimizer);
                }
                disposedValue = true;
            }
        }

        ~DFIntorBase()
        {
            Dispose(disposing: false);
        }

        public void Dispose()
        {
            Dispose(disposing: true);
            GC.SuppressFinalize(this);
        }
    }
}
