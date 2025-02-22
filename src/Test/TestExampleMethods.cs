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

        public string? BasisName { get; set; } = "sto-3g";

        [TestMethod]
        public void TestHF()
        {
            List<Atom> H2O = [
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
            List<Atom> CH3OH = [
new Atom("C", 0.664617, 0.033231, 0.000000),
new Atom("H", 1.338596, -1.873145, 0.000000),
new Atom("H", 1.338631, 0.986406, -1.650963),
new Atom("H", -1.357391, 0.033256, 0.000000),
new Atom("O", 1.565402, 1.307100, 2.206427),
new Atom("H", 0.961143, 3.017646, 2.207215),
                ];
            BasisName = "sto-3g";
            RHF rhf = new(CH3OH, BasisName);
            var totalEnergy = rhf.Run();
            Assert.AreEqual(-113.54406552303, totalEnergy, 1e-8);
            var grad = new GradRHF(rhf);
            var gradTensor = grad.Run();
            PrintForce(gradTensor, grad.Natm);
        }

        public static void PrintForce(Tensor<double> force, int natm)
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
