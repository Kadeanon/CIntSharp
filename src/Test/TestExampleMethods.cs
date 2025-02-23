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
            List<Atom> H2O = [
                new Atom("O", 0.0, 0.0, 0.0),
                new Atom("H", 0.0, 0.0, 1.0),
                new Atom("H", 0.0, 1.0, 0.0)
                ];
            BasisName = "sto-3g";
            RHF rhf = new(H2O, BasisName);
            var totalEnergy = rhf.Run();
            Assert.AreEqual(-73.45549594, totalEnergy, 1e-8);
            Console.WriteLine();
            Console.WriteLine("Calculating Gradient...");
            var grad = new GradRHF(rhf);
            var gradTensor = grad.Run();
            Console.WriteLine("Total Gradient:");
            PrintForce(gradTensor, grad.Natm);
            SequenceCloseTo(
                [0.0000000000,  2.8214032546,  2.8214032546,
                -0.0000000000,  0.2173439382, -3.0387471928,
                 0.0000000000, -3.0387471928,  0.2173439382],
                gradTensor,
                 1e-7);
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

        public static void SequenceCloseTo(
            IEnumerable<double> expected, 
            IEnumerable<double> actual,
            double tolerance)
        {
            var aEnum = expected.GetEnumerator();
            var bEnum = actual.GetEnumerator();
            while (aEnum.MoveNext() && bEnum.MoveNext())
            {
                Assert.AreEqual(aEnum.Current, bEnum.Current, tolerance);
            }
        }
    }
}
