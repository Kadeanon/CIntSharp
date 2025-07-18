using CintSharp;
using CintSharp.BasisParser;
using CintSharp.DataStructures;
using CintSharp.DataStructures.Native;
using System.Diagnostics.CodeAnalysis;
using System.Runtime.InteropServices;

namespace CintSharp.DensityFitting
{
    public class DFBuilder
    {
        public List<(Atom atom, string basisName)> Atoms { get; }

        public BasisDict BasisDict { get; }

        public Atm[]? Atms { get; private set; }

        public List<Bas>? Bases { get; private set; }

        public List<double> Envs { get; }

        public EnvBuilder PrimBuilder { get; }

        /// <summary>
        /// Create an environment builder.
        /// </summary>
        /// <param name="atoms">Atoms to build environment</param>
        /// <param name="basisParser"></param>
        public DFBuilder(EnvBuilder primBuilder, IBasisParser? basisParser = null)
        {
            PrimBuilder = primBuilder;
            Atoms = new(primBuilder.Atoms);
            Envs = primBuilder.Envs;
            BasisDict = basisParser == null ?
            primBuilder.BasisDict : new(basisParser, Envs);
        }

        /// <summary>
        /// Set basis for all atoms.
        /// </summary>
        /// <param name="basisName">The name of basis to set</param>
        /// <returns>The <see cref="EnvBuilder"/> object self</returns>
        public DFBuilder SetBasis(string basisName)
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
        public DFBuilder SetBasis(string basisName, params ReadOnlySpan<string> symbols)
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
        public DFBuilder SetBasis(string basisName, params ReadOnlySpan<int> idxs)
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
        public DFBuilder SetBasis(Func<Atom, string> basisNameSetter)
        {
            var span = CollectionsMarshal.AsSpan(Atoms);
            foreach (ref var entry in span)
            {
                entry.basisName = basisNameSetter(entry.atom);
            }
            return this;
        }

        [MemberNotNull(nameof(Atms), nameof(Bases))]
        public DFEnvs Build()
        {
            if (PrimBuilder.Atms is not null)
                throw new InvalidOperationException(
                    "The environment has been built before. " +
                    "Please create a new builder to build again.");
            var envs = PrimBuilder.Build();
            var nao = envs.NAO;
            var shellRanges = new List<Range>(Atoms.Count);
            var atomRanges = new List<Range>(Atoms.Count);
            var shellLength = new List<int>(Atoms.Count);
            Atms = envs.Atms;
            Bases = PrimBuilder.Bases;
            int current = 0;
            int ibas = 0;
            int lastOffsetByShell = 0;
            int lastOffsetByAtom = 0;
            for (int iatom = 0; iatom < Atoms.Count; iatom++)
            {
                (var atom, string basisName) = Atoms[iatom];
                var shells = BasisDict[basisName, atom.AtomNumber];
                foreach (var shell in shells)
                {
                    var shellBas = shell.CopyToBasis(iatom);
                    Bases.Add(shellBas);
                    var cgto = shellBas.CgtoSpheric();
                    shellLength.Add(cgto);
                    current += cgto;
                    shellRanges.Add(new Range(lastOffsetByShell,
                        current));
                    lastOffsetByShell = current;
                    ibas++;
                }
                atomRanges.Add(new Range(lastOffsetByAtom, current));
                lastOffsetByAtom = current;
            }

            return new DFEnvs(envs, [.. Bases], [.. Envs],
                [.. shellLength], [.. atomRanges], [.. shellRanges]);
        }
    }
}
