using CintSharp.BasisParser;
using CintSharp.Intor;
using CIntSharp.DataStructures.Native;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
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

        public IntorBase Intor(string intorName) 
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

        public Tensor<double> GetOvlp() => Intor("int1e_ovlp").Invoke();

        public Tensor<double> GetKin() => Intor("int1e_kin").Invoke();

        public Tensor<double> GetNuc() => Intor("int1e_nuc").Invoke();

        public Tensor<double> GetERI() => Intor("int2e").Invoke();

        public Tensor<double> GetHCore() => 
            Tensor.Add(GetKin().AsReadOnlyTensorSpan(), GetNuc());
    }
}
