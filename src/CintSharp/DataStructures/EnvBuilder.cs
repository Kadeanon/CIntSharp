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
        public List<(Atom atom, string basisName)> Atoms { get; }

        public BasisDict BasisDict { get; }

        public Atm[]? Atms { get; private set; }

        public List<Bas>? Bases { get; private set; }

        public List<double> Envs { get; }

        /// <summary>
        /// Create an environment builder.
        /// </summary>
        /// <param name="atoms">Atoms to build environment</param>
        /// <param name="basisParser"></param>
        public EnvBuilder(IEnumerable<Atom> atoms, IBasisParser? basisParser = null)
        {
            Atoms = atoms.Select(atom => (atom, string.Empty)).ToList();
            Envs = new(20 + atoms.Count() * 20);
            CollectionsMarshal.SetCount(Envs, 20);
            BasisDict = new(basisParser ?? new JsonBasisParser(), Envs);
        }

        /// <summary>
        /// Set basis for all atoms.
        /// </summary>
        /// <param name="basisName">The name of basis to set</param>
        /// <returns>The <see cref="EnvBuilder"/> object self</returns>
        public EnvBuilder SetBasis(string basisName)
        {
            foreach (ref var entry in CollectionsMarshal.AsSpan(Atoms))
            {
                entry.basisName = basisName;
            }
            return this;
        }

        /// <summary>
        /// Set basis for atoms with specified symbols.
        /// </summary>
        /// <param name="basisName">The name of basis to set</param>
        /// <param name="symbols">The symbols of the atoms</param>
        /// <returns>The <see cref="EnvBuilder"/> object self</returns>
        public EnvBuilder SetBasis(string basisName, params ReadOnlySpan<string> symbols)
        {
            foreach (ref var entry in CollectionsMarshal.AsSpan(Atoms))
            {
                if (symbols.Contains(entry.atom.Element.Symbol))
                {
                    entry.basisName = basisName;
                }
            }
            return this;
        }

        /// <summary>
        /// Set basis for atoms with specified indexes.
        /// </summary>
        /// <param name="basisName">The name of basis to set</param>
        /// <param name="idxs">The indexes of the atoms</param>
        /// <returns>The <see cref="EnvBuilder"/> object self</returns>
        public EnvBuilder SetBasis(string basisName, params ReadOnlySpan<int> idxs)
        {
            var span = CollectionsMarshal.AsSpan(Atoms);
            foreach (var index in idxs)
            {
                span[index].basisName = basisName;
            }
            return this;
        }

        /// <summary>
        /// Set basis with a custom setter.
        /// </summary>
        /// <param name="basisNameSetter">A custom </param>
        /// <returns></returns>
        public EnvBuilder SetBasis(Func<Atom, string> basisNameSetter)
        {
            var span = CollectionsMarshal.AsSpan(Atoms);
            foreach (ref var entry in span)
            {
                entry.basisName = basisNameSetter(entry.atom);
            }
            return this;
        }

        [MemberNotNull(nameof(Atms), nameof(Bases))]
        public CIntEnvs Build()
        {
            var shellRanges = new List<Range>(Atoms.Count);
            var atomRanges = new List<Range>(Atoms.Count);
            var shellLength = new List<int>(Atoms.Count);
            Atms = new Atm[Atoms.Count];
            Bases = new(Atoms.Count);
            int current = 0;
            int ibas = 0;
            int lastOffsetByShell = 0;
            int lastOffsetByAtom = 0;
            for(int iatom = 0; iatom < Atoms.Count; iatom++)
            {
                (Atom atom, string basisName) = Atoms[iatom];
                Shells shells = BasisDict[basisName, atom.AtomNumber];
                foreach (var shell in shells)
                {
                    Bas shellBas = new();
                    shell.CopyToBasis(ref shellBas, iatom, Envs);
                    Bases.Add(shellBas);
                    var cgto = shellBas.CgtoSpheric();
                    shellLength.Add(cgto);
                    current += cgto;
                    shellRanges.Add(new Range(lastOffsetByShell, current));
                    lastOffsetByShell = current;
                    ibas++;
                }
                atom.CopyToAtm(ref Atms[iatom], Envs);
                atomRanges.Add(new Range(lastOffsetByAtom, current));
                lastOffsetByAtom = current;
            }

            return new CIntEnvs(Atms, [.. Bases], [.. Envs], 
                [..shellLength], [..atomRanges],[..shellRanges]);
        }
    }
}
