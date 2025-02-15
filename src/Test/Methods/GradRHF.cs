using CintSharp.DataStructures;
using CintSharp.Intor;
using MathNet.Numerics.LinearAlgebra;
using NCDK.Numerics;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Test.Methods
{
    internal class GradRHF
    {
        public RHF RHF { get; set; }

        public double Energy { get; set; }

        public Tensor<double> force;

        public GradRHF(List<Atom> atoms, string basisName)
        {
            RHF = new RHF(atoms, basisName);
            Energy = RHF.Run();
            force = Tensor.CreateUninitialized<double>([Natm * 3]);
        }

        public GradRHF(RHF rhf)
        {
            RHF = rhf;
            Energy = rhf.TotalEnergy;
            force = Tensor.CreateUninitialized<double>([Natm * 3]);
        }

        public Tensor<double> Run()
        {
            CalculateForce();
            return force;
        }


        public Tensor<double> overlapDerivative;
        public Tensor<double> hamiltonianDerivative;
        public Tensor<double> eriDerivative;
        public Tensor<double> elecEnergyDerivative;
        public Tensor<double> nuclearRepulsionDerivative;
        public Tensor<double> overlapDerivativeMO;

        public int Nao => RHF.Envs.NAO;

        public int Natm => RHF.Envs.Natm;

        public void CalculateForce()
        {
            #region prepare
            CIntEnvs mol = RHF.Envs;
            var fock = Matrix2Tensor(RHF.FockMatrix);
            var coeff = Matrix2Tensor(RHF.C);
            var density = Matrix2Tensor(RHF.P);//GradientUtils.AO2MO(coeff, Core.P);
            var nocc = RHF.Atoms.Sum(atm => atm.AtomNumber) / 2;
            var so = new Range(0, nocc);
            var ac = RHF.Atoms.Select(atm => atm.position).ToArray();
            var nao = Nao;
            var length = Natm * 3;

            var int1e_ipovlp = mol.InvokeIntor("int1e_ipovlp");
            var int1e_ipkin = mol.InvokeIntor("int1e_ipkin");
            var int1e_ipnucl = mol.InvokeIntor("int1e_ipnuc");
            var z_a = RHF.Atoms.Select(atm => atm.AtomNumber).ToArray();
            #endregion
            #region overlapDerivativ and hamiltonianDerivative
            //derivatives in AO basis
            overlapDerivative = Tensor.Create<double>([length, nao, nao]);
            hamiltonianDerivative = Tensor.Create<double>([length, nao, nao]);
            for (var iatm = 0; iatm < Natm; iatm++)
            {
                var p0 = mol.OffsetsByAtoms[iatm];
                var p1 = mol.OffsetsByAtoms[iatm + 1];
                var sa = Enumerable.Range(p0, p1 - p0).ToList();
                Tensor<double> int1e_iprinv;
                using (mol.RinvAt(iatm)) 
                {
                    int1e_iprinv = mol.InvokeIntor("int1e_iprinv");
                }
                for (int idim = 0; idim < 3; idim++)
                {
                    int totalDim = iatm * 3 + idim;
                    for (var k = 0; k < nao; k++)
                    {
                        foreach (var j in sa)
                        {
                            double value = int1e_ipovlp[j, k, idim];
                            overlapDerivative[totalDim, j, k] -= value;
                            overlapDerivative[totalDim, k, j] -= value;
                            value = int1e_ipkin[j, k, idim] + int1e_ipnucl[j, k, idim];
                            hamiltonianDerivative[totalDim, j, k] -= value;
                            hamiltonianDerivative[totalDim, k, j] -= value;
                        }
                        for (int j = 0; j < nao; j++)
                        {
                            double value = z_a[iatm] * int1e_iprinv[j, k, idim];
                            hamiltonianDerivative[totalDim, j, k] -= value;
                            hamiltonianDerivative[totalDim, k, j] -= value;
                        }
                    }
                }
            }
            #endregion
            #region eriDerivative
            var int2e_ip1 = mol.InvokeIntor("int2e_ip1");
            //var sum_int2e_ip1 = int2e_ip1.square().sum().ToDouble();
            var eri1_ao = Tensor.Create<double>([3 * Natm, Nao, Nao, Nao, Nao]);
            for (var iatm = 0; iatm < Natm; iatm++)
            {
                int p0 = mol.OffsetsByAtoms[iatm];
                int p1 = mol.OffsetsByAtoms[iatm + 1];
                for (int idim = 0; idim < 3; idim++)
                {
                    int totalDim = iatm * 3 + idim;
                    for (var i = 0; i < nao; i++)
                    {
                        for (var j = 0; j < nao; j++)
                        {
                            for (var k = 0; k < nao; k++)
                            {
                                for (var s = p0; s < p1; s++)
                                {
                                    eri1_ao[totalDim, s, i, j, k] -= int2e_ip1[s, i, j, k, idim];
                                    eri1_ao[totalDim, i, s, j, k] -= int2e_ip1[s, i, j, k, idim];
                                    eri1_ao[totalDim, i, j, s, k] -= int2e_ip1[s, k, i, j, idim];
                                    eri1_ao[totalDim, i, j, k, s] -= int2e_ip1[s, k, i, j, idim];
                                }
                            }
                        }
                    }
                }
            }
            //sum_eri = eri1_ao.square().sum().ToDouble();//83.409308874986692
            //var sum_eri_1 = eri1_ao.abs().sum().ToDouble();//1594.7648127209497
            eriDerivative = eri1_ao;
            #endregion
            #region electronicEnergyDerivative
            overlapDerivativeMO = Tensor.Create<double>([length, Nao, Nao]);
            var fock_mo = AO2MO(coeff, fock);
            for (int itotalDim = 0; itotalDim < length; itotalDim++)
            {
                Span<NRange> ranges = [ new(itotalDim, itotalDim + 1), NRange.All, NRange.All];
                var span = overlapDerivative.Slice(ranges);
                span.SetSlice(AO2MO(coeff, span));
            }

            elecEnergyDerivative = Tensor.Create<double>([length]);
            var sub1 = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(length);
            var sub2 = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(length);
            var sub3 = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(length);
            var sub4 = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(length);
            var so_enum = Enumerable.Range(0, nocc).ToArray();
            for (var iatm = 0; iatm < Natm; iatm++)
            {
                int p0 = mol.OffsetsByAtoms[iatm];
                int p1 = mol.OffsetsByAtoms[iatm + 1];
                var slice = Enumerable.Range(p0, p1 - p1).ToArray();
                //sub1 += torch.einsum("tij,ij -> t", hamiltonianDerivative, density);
                for (int t = 0; t < length; t++) 
                {
                    double it = 0;
                    for (int i = 0; i < nao; i++)
                    {
                        for (int j = 0; j < nao; j++)
                        {
                            it += hamiltonianDerivative[t, i, j] * density[i, j];
                        }
                    }
                    sub1[t] += it;
                }
            }
            //sub2 = torch.einsum("duvkl,uv,kl -> d", eriDerivative, density, density);
            var dkl = Tensor.Create<double>([length, Nao, Nao]);
            var duvkl = eriDerivative;
            var uv = density;
            for (int d = 0; d < length; d++)
            {
                for (int k = 0; k < Nao; k++)
                {
                    for (int l = 0; l < Nao; l++)
                    {
                        double idkl = 0;
                        for (int u = 0; u < Nao; u++)
                        {
                            for (int v = 0; v < Nao; v++)
                            {
                                idkl += duvkl[d, u, v, k, l] * uv[u, v];
                            }
                        }
                        dkl[d, k, l] = idkl;
                    }
                }
            }
            var kl = density;
            for (int d = 0; d < length; d++)
            {
                double id = 0;
                for (int k = 0; k < Nao; k++)
                {
                    for (int l = 0; l < Nao; l++)
                    {
                        id += kl[k, l] * dkl[d, k, l];
                    }
                }
                sub2[d] = id;
            }
            //sub3 = -torch.einsum("dukvl,uv,kl -> d", eriDerivative, density, density);
            var dukvl = eriDerivative;
            for (int d = 0; d < length; d++)
            {
                for(int k = 0; k < Nao; k++) 
                {
                    for(int l = 0; l < Nao; l++)
                    {
                        double idkl = 0;
                        for (int u = 0; u < Nao; u++)
                        {
                            for (int v = 0; v < Nao; v++)
                            {
                                idkl += dukvl[d, u, k, v, l] * uv[u, v];
                            }
                        }
                        dkl[d, k, l] = idkl;
                    }
                }
            }
            for (int d = 0; d < length; d++)
            {
                double id = 0;
                for (int k = 0; k < Nao; k++)
                {
                    for (int l = 0; l < Nao; l++)
                    {
                        id += kl[k, l] * dkl[d, k, l];
                    }
                }
                sub3[d] = id;
            }
            //sub4 -= torch.einsum("dij, ij -> d", overlapDerivativeMO, fock_mo);
            var dij = overlapDerivativeMO;
            var ij = fock_mo;
            for (int d = 0; d < length; d++)
            {
                double id = 0;
                for (int i = 0; i < Nao; i++)
                {
                    for (int j = 0; j < Nao; j++)
                    {
                        id += dij[d, i, j] * ij[i, j];
                    }
                }
                sub4[d] = id;
            }
            elecEnergyDerivative = Vector2Tensor(sub1 + sub2 / 2 - sub3 / 4 - sub4 * 2);
            nuclearRepulsionDerivative = Tensor.Create<double>([length]);
            #endregion
            #region nuclear repulsion derivative
            for (var iatm = 0; iatm < Natm; iatm++)
            {
                var z_i = z_a[iatm];
                var coor_i = ac[iatm];
                for (int idim = 0; idim < 3; idim++)
                {
                    int totalDim = iatm * 3 + idim;
                    double value = 0;
                    for (var jatm = 0; jatm < Natm; jatm++)
                    {
                        if (jatm == iatm)
                        {
                            continue;
                        }

                        var z_j = z_a[jatm];
                        var coor_j = ac[jatm];
                        var coor = coor_i - coor_j;
                        var coorSpan = MemoryMarshal.CreateSpan(ref Unsafe.As<Vector3, double>(ref coor), 3);
                        var z_ij = z_i * z_j;
                        value -= z_ij * coorSpan[idim] / Math.Pow(coor.Length(), 3);
                    }

                    nuclearRepulsionDerivative[totalDim] = value;
                }
            }
            #endregion
            force = Tensor.Add(elecEnergyDerivative.AsReadOnlyTensorSpan(), nuclearRepulsionDerivative);
        }

        private static Tensor<double> AO2MO(ReadOnlyTensorSpan<double> moCoeffs, ReadOnlyTensorSpan<double> getter)
        {
            var dims = getter.Lengths;
            if (moCoeffs.Rank != 2) 
            {
                throw new ArgumentException("The input matrix must be 2D.");
            }
            nint length = moCoeffs.Lengths[0];
            if (length != moCoeffs.Lengths[1])
            {
                throw new ArgumentException("The input matrix must be square.");
            }
            foreach (var dim in dims)
            {
                if (dim != length)
                {
                    throw new ArgumentException("The input matrix must have the same size as the moCoeffs.");
                }
            }

            Tensor<double> temp1 = Tensor.Create<double>(dims);
            Tensor<double> temp2 = Tensor.Create<double>(dims);

            switch (dims.Length) 
            {
                case 2:
                    {
                        var uv = getter;
                        var up = moCoeffs;
                        var pv = temp1;
                        for (int p = 0; p < length; p++)
                        {
                            for (int v = 0; v < length; v++)
                            {
                                var ipv = 0.0;
                                for (int u = 0; u < length; u++)
                                {
                                    ipv += uv[u, v] * up[u, p];
                                }
                                pv[p, v] = ipv;
                            }
                        }
                        var vq = moCoeffs;
                        var pq = temp2;
                        for (int p = 0; p < length; p++)
                        {
                            for (int q = 0; q < length; q++)
                            {
                                var ipq = 0.0;
                                for (int v = 0; v < length; v++)
                                {
                                    ipq += pv[p, v] * vq[v, q];
                                }
                                pq[p, q] = ipq;
                            }
                        }
                        return pq;
                    }
                case 4:
                    {
                        var uvkl = getter;
                        var up = moCoeffs;
                        var pvkl = temp1;
                        for (int p = 0; p < length; p++)
                        {
                            for (int v = 0; v < length; v++)
                            {
                                for (int k = 0; k < length; k++)
                                {
                                    for (int l = 0; l < length; l++)
                                    {
                                        var ipvkl = 0.0;
                                        for (int u = 0; u < length; u++)
                                        {
                                            ipvkl += uvkl[u, v, k, l] * up[u, p];
                                        }
                                        pvkl[p, v, k, l] = ipvkl;
                                    }
                                }
                            }
                        }
                        var vq = moCoeffs;
                        var pqkl = temp2;
                        for (int p = 0; p < length; p++)
                        {
                            for (int q = 0; q < length; q++)
                            {
                                for (int k = 0; k < length; k++)
                                {
                                    for (int l = 0; l < length; l++)
                                    {
                                        var ipqkl = 0.0;
                                        for (int v = 0; v < length; v++)
                                        {
                                            ipqkl += pvkl[p, v, k, l] * vq[v, q];
                                        }
                                        pqkl[p, q, k, l] = ipqkl;
                                    }
                                }
                            }
                        }
                        var kr = moCoeffs;
                        var pqrl = temp1;
                        for (int p = 0; p < length; p++)
                        {
                            for (int q = 0; q < length; q++)
                            {
                                for (int r = 0; r < length; r++)
                                {
                                    for (int l = 0; l < length; l++)
                                    {
                                        var ipqrl = 0.0;
                                        for (int k = 0; k < length; k++)
                                        {
                                            ipqrl += pqkl[p, q, k, l] * kr[k, r];
                                        }
                                        pqrl[p, q, r, l] = ipqrl;
                                    }
                                }
                            }
                        }
                        var ls = moCoeffs;
                        var pqrs = temp2;
                        for (int p = 0; p < length; p++)
                        {
                            for (int q = 0; q < length; q++)
                            {
                                for (int r = 0; r < length; r++)
                                {
                                    for (int s = 0; s < length; s++)
                                    {
                                        var ipqrs = 0.0;
                                        for (int l = 0; l < length; l++)
                                        {
                                            ipqrs += pqrl[p, q, r, l] * ls[l, s];
                                        }
                                        pqrs[p, q, r, s] = ipqrs;
                                    }
                                }
                            }
                        }
                        return pqrs;
                    }
                default:
                    throw new ArgumentException("The input matrix must be 2D or 4D.");
            }
        }

        private Tensor<double> Matrix2Tensor(Matrix<double> matrix)
        {
            Tensor<double> tensor = Tensor.CreateUninitialized<double>([matrix.RowCount, matrix.ColumnCount]);
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    tensor[i, j] = matrix[i, j];
                }
            }
            return tensor;
        }

        private Tensor<double> Vector2Tensor(MathNet.Numerics.LinearAlgebra.Vector<double> matrix)
        {
            Tensor<double> tensor = Tensor.CreateUninitialized<double>([matrix.Count]);
            for (int i = 0; i < matrix.Count; i++)
            {
                tensor[i] = matrix[i];
            }
            return tensor;
        }
    }
}
