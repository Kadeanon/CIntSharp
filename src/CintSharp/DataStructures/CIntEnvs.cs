using CintSharp.BasisParser;
using CintSharp.DataStructures.Native;
using CintSharp.Intor;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DataStructures
{
    public class CIntEnvs(Atm[] atms, Bas[] bases, double[] envs, int[] offsetsByAtoms, int[] offsetsByShells)
    {

        public Atm[] Atms { get; } = atms;

        public Bas[] Bases { get; } = bases;

        public double[] Envs { get; } = envs;

        public int[] OffsetsByAtoms { get; } = offsetsByAtoms;

        public int[] OffsetsByShells { get; } = offsetsByShells;

        public ref EnvHeader EnvHeader => ref EnvHeader.FromSpan(Envs);

        public int Natm => Atms.Length;

        public int Nbas => Bases.Length;

        public int NAO => OffsetsByAtoms[^1];

        public static CIntEnvs Create(List<Atom> atoms, Func<Atom, string> basisNameSetter, IBasisParser parser = null) 
            => new EnvBuilder(atoms, parser).SetBasis(basisNameSetter).Build();

        public IntorBase CreateIntor(string intorName) 
        {
            if (intorName.StartsWith("int1e"))
            {
                return new Intor1e(this, intorName);
            }
            else if (intorName.StartsWith("int2e"))
            {
                return new Intor2e(this, intorName);
            }
            else 
            {
                throw new NotImplementedException();
            }
        }

        public Vector3 AtomCoord(int atomIndex) 
        {
            return Unsafe.As<double, Vector3>(ref Envs[Atms[atomIndex].pointerForCoords]);
        }

        public HeaderRinvScope RinvAt(int atomIndex)
        {
            return new HeaderRinvScope(this, atomIndex);
        }

        public HeaderRinvScope RinvAt(Vector3 rinvOrigin)
        {
            return new HeaderRinvScope(this, rinvOrigin);
        }
    }

    public readonly struct HeaderRinvScope : IDisposable 
    {
        public CIntEnvs Envs { get; }

        public readonly ref Vector3 Rinv => ref Envs.EnvHeader.RinvOrigin;

        public Vector3 RinvOrigin { get; }

        public HeaderRinvScope(CIntEnvs envs, Vector3 rinvOrigin)
        {
            Envs = envs;
            ref Vector3 rinv = ref Envs.EnvHeader.RinvOrigin;
            RinvOrigin = rinv;
            rinv = rinvOrigin;
        }

        public HeaderRinvScope(CIntEnvs envs, int atomIndex) : this(envs, envs.AtomCoord(atomIndex))
        {
        }

        public void Dispose()
        {
            Envs.EnvHeader.RinvOrigin = RinvOrigin;
        }
    }
}
