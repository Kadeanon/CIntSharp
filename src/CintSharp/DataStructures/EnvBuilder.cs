using CintSharp.BasisParser;
using CintSharp.DataStructures;
using CintSharp.Native.Libcint;
using CintSharp.DataStructures.Native;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DataStructures
{
    public class EnvBuilder
    {
        public List<Atom> Atoms { get; }

        public BasisDict BasisDict { get; }

        public Atm[] Atms { get; }

        public List<Bas> Bases { get; }

        public List<double> Envs { get; }


        public EnvBuilder(List<Atom> atoms, IBasisParser basisParser = null)
        {
            Atoms = atoms;
            Atms = new Atm[atoms.Count];
            Bases = new(atoms.Count);
            Envs = new(20 + atoms.Count * 20);
            CollectionsMarshal.SetCount(Envs, 20);
            BasisDict = new BasisDict(basisParser ?? new JsonBasisParser(), Envs);
        }

        public EnvBuilder SetBasis(string basisName)
        {
            Span<Atom> span = CollectionsMarshal.AsSpan(Atoms);
            foreach (ref Atom atom in span)
            {
                atom.BasisName = basisName;
            }
            return this;
        }

        public EnvBuilder SetBasis(string basisName, params ReadOnlySpan<string> symbols)
        {
            Span<Atom> span = CollectionsMarshal.AsSpan(Atoms);
            foreach (ref Atom atom in span)
            {
                if (symbols.Contains(atom.Element.Symbol))
                {
                    atom.BasisName = basisName;
                }
            }
            return this;
        }

        public EnvBuilder SetBasis(string basisName, params ReadOnlySpan<int> idxs)
        {
            Span<Atom> span = CollectionsMarshal.AsSpan(Atoms);
            foreach (var index in idxs)
            {
                span[index].BasisName = basisName;
            }
            return this;
        }

        public EnvBuilder SetBasis(Func<Atom, string> basisNameSetter)
        {
            Span<Atom> span = CollectionsMarshal.AsSpan(Atoms);
            foreach (ref Atom atom in span)
            {
                atom.BasisName = basisNameSetter(atom);
            }
            return this;
        }

        public CIntEnvs Build()
        {
            List<int> offsetsByAtom = [];
            List<int> offsetsByShell = [];
            int current = 0;
            int ibas = 0;
            int atomIdx = 0;
            foreach (var atom in Atoms)
            {
                Shells shells = BasisDict[atom];
                offsetsByAtom.Add(current);
                foreach (var shell in shells)
                {
                    Bas shellBas = new();
                    offsetsByShell.Add(current);
                    shell.CopyToBasis(ref shellBas, atomIdx, Envs);
                    Bases.Add(shellBas);
                    current += shellBas.CgtoSpheric();
                    ibas++;
                }
                atom.CopyToAtm(ref Atms[atomIdx], Envs);
                atomIdx++;
            }

            offsetsByAtom.Add(current);
            offsetsByShell.Add(current);

            return new CIntEnvs(Atms, Bases.ToArray(), Envs.ToArray(), 
                [..offsetsByAtom], [..offsetsByShell]);
        }
    }
}
