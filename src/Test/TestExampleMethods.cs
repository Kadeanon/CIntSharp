using CintSharp.DataStructures;
using CintSharp.Test.Methods;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Numerics.Tensors;
using System.Text;
using System.Threading.Tasks;

using Vector = MathNet.Numerics.LinearAlgebra.Vector<double>;

namespace CintSharp.Test
{
    [TestClass]
    public class TestExampleMethods
    {
        public List<Atom>? H2O { get; set; }

        public string? BasisName { get; set; }

        [TestMethod]
        public void TestHF()
        {
            H2O = [
                new Atom("O", 0.0, 0.0, 0.0),
                new Atom("H", 0.0, 0.0, 1.0),
                new Atom("H", 0.0, 1.0, 0.0)
                ];
            BasisName = "sto-3g";
            RHF rhf = new(H2O, BasisName);
            var totalEnergy = rhf.Run();
            Assert.AreEqual(-73.45549594, totalEnergy, 1e-8);
        }

        //[Ignore("This example code is not implemented yet.")]
        [TestMethod]
        public void TestGrad()
        {
            H2O = [
                new Atom("O", 0.0, 0.0, 0.0),
                new Atom("H", 0.0, 0.0, 1.0),
                new Atom("H", 0.0, 1.0, 0.0)
                ];
            BasisName = "sto-3g";
            RHF rhf = new(H2O, BasisName);
            var totalEnergy = rhf.Run();
            Assert.AreEqual(-73.45549594, totalEnergy, 1e-8);
            var grad = new GradRHF(rhf);
            var gradTensor = grad.Run();
            PrintForce(gradTensor, grad.Natm);
        }

        public void PrintForce(Tensor<double> force, int natm)
        {
            StringBuilder sb = new();
            for (int i = 0; i < natm; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sb.Append($"{force[i * 3 + j]:F6} ");
                }
                sb.AppendLine();
            }
            Console.WriteLine(sb.ToString());
        }
    }
}
