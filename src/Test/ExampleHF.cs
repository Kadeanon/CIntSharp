using CintSharp.DataStructures;
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
    public class ExampleHF
    {
        public List<Atom> Atoms { get; set; }

        public CIntEnvs H2OEnvs { get; set; }

        public Matrix<double> Ovlp { get; set; }

        public Matrix<double> HCore { get; set; }

        public Tensor<double> ERI { get; set; }

        public Matrix<double> X { get; set; }

        public Matrix<double> XH { get; set; }

        public Matrix<double> P { get; set; }

        public Matrix<double> G{get; set;}

        public Matrix<double> C { get; set; }

        public Matrix<double> FockMatrix{get; set;}

        public double[] es { get; set; }

        public const int MaxCycle = 1024;
        public const double RequiredDelEnergy = 1e-8;
        private const double RequiredMaxDelF = 1e-8;
        private const double RequiredAverDelF = 1e-8;

        public double ElectronicEnergy
        {
            get
            {
                int nao = H2OEnvs.NAO;
                double d = 0;
                for (int i = 0; i < nao; i++)
                {
                    for (int j = 0; j < nao; j++)
                    {
                        d += P[j, i] * (HCore[i, j] + FockMatrix[i, j]);
                    }
                }
                return d / 2;
            }
        }

        public double TotalEnergy
        {
            get
            {
                double result = ElectronicEnergy;
                for (int i = 0; i < Atoms.Count; i++)
                {
                    for (int j = i + 1; j < Atoms.Count(); j++)
                    {
                        var a = Atoms[i];
                        var b = Atoms[j];
                        result += a.AtomNumber * b.AtomNumber / (a.position - b.position).Length();
                    }
                }
                return result;
            }
        }


        [TestMethod]
        public void RunHF() 
        {
            InitExample();
            Conver();
            Assert.AreEqual(-73.45549594, TotalEnergy, 1e-8);
        }

        public void InitExample()
        {
            CreateH2OEnvs();
            int nao = H2OEnvs.NAO;
            Ovlp = Tensor2Matrix(H2OEnvs.GetOvlp());
            HCore = Tensor2Matrix(H2OEnvs.GetHCore());
            ERI = H2OEnvs.GetERI();
            X = ReverseSqrt(Ovlp);
            XH = X.Transpose();
            P = Matrix<double>.Build.Dense(nao, nao);
            GeneralGMatrix();
            FockMatrix = HCore + G;
        }

        public void CreateH2OEnvs() 
        {
            Atoms = [
                new Atom("O", 0.0, 0.0, 0.0),
                new Atom("H", 0.0, 0.0, 1.0),
                new Atom("H", 0.0, 1.0, 0.0)
                ];
            H2OEnvs = CIntEnvs.Create(Atoms, atm => "sto-3g");
        }

        public static Matrix<double> Tensor2Matrix(Tensor<double> tensor)
        {
            if (tensor.Rank != 2)
            {
                throw new ArgumentException("The tensor must be a matrix.");
            }
            Matrix<double> matrix = Matrix<double>.Build.Dense((int)tensor.Lengths[0], (int)tensor.Lengths[1]);
            for (int i = 0; i < tensor.Lengths[0]; i++)
            {
                for (int j = 0; j < tensor.Lengths[1]; j++)
                {
                    matrix[i, j] = tensor[i, j];
                }
            }
            return matrix;
        }

        public static Matrix<double> ReverseSqrt(Matrix<double> Data)
        {
            var copy = Data;
            var evd = copy.Evd(Symmetricity.Hermitian);
            var eigenValues = evd.EigenValues;
            var U = evd.EigenVectors;
            var sss = U.Transpose() * Data * U;
            var s = Matrix<double>.Build.Diagonal(Data.RowCount, Data.ColumnCount, x => 1 / Math.Sqrt(eigenValues[x].Real));
            var result = U * s * U.Transpose();
            return result;
        }

        public void GeneralGMatrix()
        {
            int basisLength = P.RowCount;
            G = Matrix<double>.Build.SameAs(P);
            for (int u = 0; u < basisLength; u++)
            {
                for (int v = 0; v < basisLength; v++)
                {
                    double value = 0;
                    for (int sigma = 0; sigma < basisLength; sigma++)
                    {
                        for (int lambda = 0; lambda < basisLength; lambda++)
                        {
                            double item = P[sigma, lambda] * (2 * ERI[u, v, sigma, lambda] - ERI[u, sigma, v, lambda]) / 2;
                            value += item;
                        }
                    }
                    G[u, v] = value;
                }
            }
        }

        public void GeneralPMatrix()
        {
            int occ = Atoms.Sum(atom => atom.AtomNumber) / 2;
            P = Matrix<double>.Build.SameAs(C);

            for (int i = 0; i < C.ColumnCount; i++)
            {
                for (int j = 0; j < C.ColumnCount; j++)
                {
                    for (int index = 0; index < occ; index++)
                    {
                        P[i, j] += 2 * C[i, index] * C[j, index];
                    }
                }
            }
        }

        public void Conver()
        {
            int cycle = 0;
            bool isConver;
            do
            {
                cycle++;
                Console.WriteLine($"Cycle {cycle}:");
                double oldEnergy = ElectronicEnergy;
                var oldFockMatrix = FockMatrix;
                ConverInner(cycle);
                double delEnergy = oldEnergy - ElectronicEnergy;
                var delFock = (FockMatrix - oldFockMatrix).PointwiseAbs();
                isConver = ConjConver(delEnergy, delFock);
                Console.WriteLine();
            } while (cycle <= MaxCycle && !isConver);
            Console.WriteLine();
            if (!isConver)
            {
                Console.WriteLine($"Not Converged after Cycle {MaxCycle}!");
            }
            else
            {
                Console.WriteLine($"Converged in Cycle {cycle}! Totol Energy= {TotalEnergy}, Electronic Energy= {ElectronicEnergy}.");
            }

        }

        public void ConverInner(int cycle)
        {
            var F_ = XH * FockMatrix * X;
            var evd = F_.Evd(Symmetricity.Symmetric);
            var C_ = evd.EigenVectors;
            es = evd.EigenValues.Real().ToArray();
            C = X * C_;
            SortEvd();
            GeneralPMatrix();
            GeneralGMatrix();
            FockMatrix = HCore + G;
            Console.WriteLine($"Energy: {TotalEnergy}");
        }

        public void SortEvd()
        {
            Dictionary<double, Vector> evd = new();
            for (int i = 0; i < es.Length; i++)
            {
                evd.Add(es[i], C.Column(i));
            }
            var kvps = evd.ToList();
            kvps.Sort((x1, x2) => x1.Key.CompareTo(x2.Key));
            es = kvps.Select(x => x.Key).ToArray();
            for (int i = 0; i < es.Length; i++)
            {
                C.SetColumn(i, kvps[i].Value);
            }
            OutputEnergys();
        }

        public void OutputEnergys()
        {
            int occ = Atoms.Sum(atom => atom.AtomNumber) / 2;
            int j = 0;
            Console.WriteLine("elec energy: ");
            StringBuilder occStr = new("Occ: [");
            for (; j < occ; j++)
            {
                occStr.Append($"{es[j]:F4}");
                if (j < occ - 1) occStr.Append(", ");
            }
            occStr.Append(']');
            Console.WriteLine(occStr.ToString());
            StringBuilder virtStr = new("Virt: [");
            for (; j < es.Length; j++)
            {
                virtStr.Append($"{es[j]:F4}");
                if (j < es.Length - 1) virtStr.Append(", ");
            }
            virtStr.Append(']');
            Console.WriteLine(virtStr.ToString());
        }

        public bool ConjConver(double delEnergy, Matrix<double> delFock)
        {

            string Conj(double result, double need)
            {
                string symbol = result > need ? ">" : "<";
                if (Math.Abs(result) > 1e-4)
                {
                    return $"{result:F4}({symbol}{need})";
                }
                return $"{result:E2}({symbol}{need})";
            }
            double dropEnergy = delEnergy;
            double maxDelF = 0;
            double totalSqDelF = 0;
            for (int i = 0; i < delFock.ColumnCount; i++)
            {
                for (int j = 0; j < delFock.ColumnCount; j++)
                {
                    double d = delFock[i, j];
                    totalSqDelF += d * d;
                    if (d > maxDelF)
                    {
                        maxDelF = d;
                    }
                }
            }
            double averDelF = totalSqDelF / (P.ColumnCount * P.ColumnCount);
            bool isConver = Math.Abs(dropEnergy) < RequiredDelEnergy && maxDelF < RequiredMaxDelF && averDelF < RequiredAverDelF;

            Console.WriteLine($"Energy Drop: {Conj(dropEnergy, RequiredDelEnergy)}, Max DelF: {Conj(maxDelF, RequiredMaxDelF)}, RMS DelF: {Conj(averDelF, RequiredAverDelF)}");
            if (isConver)
            {
                Console.WriteLine("Converged!");
            }
            else
            {
                Console.WriteLine("Not Converged!");
            }
            return isConver;
        }

    }
}
