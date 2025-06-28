using CintSharp.BasisParser;
using CintSharp.DataStructures.Native;
using CintSharp.Intor;
using System;
using System.Collections.Generic;
using System.IO.Pipes;
using System.Linq;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DataStructures
{
    public class CIntEnvs
    {

        public Atm[] Atms { get; }

        public Bas[] Bases { get; }

        public double[] Envs { get; }

        public Range[] RangesByShells { get; }

        public Range[] RangesByAtoms { get; }

        public int[] ShellLengths { get; }

        public ref EnvHeader EnvHeader => ref EnvHeader.FromSpan(Envs);

        public int Natm => Atms.Length;

        public int Nbas => Bases.Length;

        public int NAO { get; }

        internal CIntEnvs(Atm[] atms, Bas[] bases, double[] envs, int[] lengthByShells, Range[] rangesByAtoms, Range[] rangesByShells)
        {
            Atms = atms;
            Bases = bases;
            Envs = envs;
            RangesByShells = rangesByShells;
            RangesByAtoms = rangesByAtoms;
            ShellLengths = lengthByShells;
            NAO = lengthByShells.Sum();
        }

        public static CIntEnvs Create(IEnumerable<Atom> atoms, Func<Atom, string> basisNameSetter, IBasisParser? parser = null) 
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

        public IEnumerable<int> EnumerateByAtom(int iatm) 
        {
            var range = RangesByAtoms[iatm];
            (int start, int length) = range.GetOffsetAndLength(NAO);
            return Enumerable.Range(start, length);
        }

        public IEnumerable<int> EnumerateByShell(int ish)
        {
            var range = RangesByShells[ish];
            (int start, int length) = range.GetOffsetAndLength(NAO);
            return Enumerable.Range(start, length);
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
