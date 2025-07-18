using CintSharp.BasisParser;
using CintSharp.DataStructures;
using CintSharp.DataStructures.Native;
using CintSharp.DensityFitting.Intor;

namespace CintSharp.DensityFitting
{
    public class DFEnvs : IEnvorinment
    {
        public Atm[] Atms { get; }

        public Bas[] Bases { get; }

        public double[] Envs { get; }

        public Range[] RangesByShells { get; }

        public Range[] RangesByAtoms { get; }

        public int[] ShellLengths { get; }

        public ref EnvHeader EnvHeader => ref EnvHeader.FromSpan(Envs);

        public int Natm => Atms.Length;

        public int Nbas { get; }

        public int NAO { get; }

        public Range[] RangesByShellsAux { get; }

        public Range[] RangesByAtomsAux { get; }

        public int[] ShellLengthsAux { get; }

        public int NAOAux { get; }

        internal DFEnvs(CIntEnvs primEnvs, Bas[] bases, double[] envs,
            int[] lengthByShellsAux, Range[] rangesByAtomsAux, Range[] rangesByShellsAux)
        {
            Atms = primEnvs.Atms;
            Nbas = bases.Length;
            RangesByShells = primEnvs.RangesByShells;
            RangesByAtoms = primEnvs.RangesByAtoms;
            ShellLengths = primEnvs.ShellLengths;
            NAO = primEnvs.NAO;
            Bases = bases;
            Envs = envs;
            NAOAux = lengthByShellsAux.Sum();
            RangesByAtomsAux = rangesByAtomsAux;
            ShellLengthsAux = lengthByShellsAux;
            RangesByShellsAux = rangesByShellsAux;
        }

        public static DFEnvs Create(
            IEnumerable<Atom> atoms, Func<Atom, string> basisNameSetter,
            Func<Atom, string> basisNameSetterAux, IBasisParser? parser = null)
        {
            var primBuilder = new EnvBuilder(atoms, parser).
                SetBasis(basisNameSetter);
            var dfBuilder = new DFBuilder(primBuilder)
                .SetBasis(basisNameSetterAux);
            return dfBuilder.Build();
        }


        public DFIntorBase CreateIntor(string intorName)
        {
            if (intorName == "int2c2e")
            {
                return new DFIntor2c2e(this);
            }
            else if (intorName == "int3c2e")
            {
                return new DFIntor3c2e(this);
            }
            else if(intorName.StartsWith("int1e"))
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
}
