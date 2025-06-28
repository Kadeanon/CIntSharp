using CintSharp.DataStructures;
using CintSharp.Intor;
using SimpleHelpers;
using SimpleHelpers.MultiAlg;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Examples;

internal class GradRHF
{
    public RHF RHF { get; set; }

    public double Energy { get; set; }

    public NDArray force;

    public GradRHF(List<Atom> atoms, string basisName)
    {
        RHF = new RHF(atoms, basisName);
        Energy = RHF.Run();
        force = NDArray.CreateUninitialized([Natm * 3]);
    }

    public GradRHF(RHF rhf)
    {
        RHF = rhf;
        Energy = rhf.TotalEnergy;
        force = NDArray.CreateUninitialized([Natm * 3]);
    }

    public NDArray Run()
    {
        CalculateForce();
        return force;
    }


    public NDArray overlapDerivative;
    public NDArray hamiltonianDerivative;
    public NDArray eriDerivative;
    public NDArray elecEnergyDerivative;
    public NDArray nuclearRepulsionDerivative;
    public NDArray overlapDerivativeMO;

    public int Nao => RHF.Envs.NAO;

    public int Natm => RHF.Envs.Natm;

    public void CalculateForce()
    {
        #region prepare
        var mol = RHF.Envs;
        var fock = RHF.Fock.AsNDArray();
        var coeff = RHF.C!.AsNDArray();
        var density = RHF.P.AsNDArray();
        var nocc = RHF.Atoms.Sum(atm => atm.AtomNumber) / 2;
        var so = new Range(0, nocc);
        var ac = RHF.Atoms.Select(atm => atm.position).ToArray();
        var nao = Nao;
        var length = Natm * 3;

        var int1e_ipovlp = mol.InvokeIntor("int1e_ipovlp");
        var int1e_ipkin = mol.InvokeIntor("int1e_ipkin");
        var int1e_ipnucl = mol.InvokeIntor("int1e_ipnuc");
        var z_a = RHF.Atoms.Select(atm => atm.AtomNumber).ToArray();
        double d49 = int1e_ipnucl[0, 3, 2];
        double d20 = int1e_ipnucl[0, 2, 3];
        double d05 = int1e_ipnucl[1, 6, 4];
        double d54 = int1e_ipnucl[1, 5, 4];

        //int1e_ipnucl[0, 2, 3] = d49;
        //int1e_ipnucl[0, 2, 4] = d49;//
        //int1e_ipnucl[0, 3, 2] = d20;//
        //int1e_ipnucl[0, 4, 2] = d20;//
        //int1e_ipnucl[1, 3, 4] = d49;//
        //int1e_ipnucl[1, 4, 3] = d20;//
        //int1e_ipnucl[1, 5, 4] = d05;//
        //int1e_ipnucl[1, 6, 4] = d54;//
        //int1e_ipnucl[2, 3, 4] = d20;//
        //int1e_ipnucl[2, 4, 3] = d49;//
        //int1e_ipnucl[2, 5, 3] = d54;//
        //int1e_ipnucl[2, 6, 3] = d05;//
        var int1e_core = int1e_ipkin + int1e_ipnucl;


        #endregion

        #region overlapDerivativ and hamiltonianDerivative
        //derivatives in AO basis
        overlapDerivative = NDArray.Create([length, nao, nao]);
        hamiltonianDerivative = NDArray.Create([length, nao, nao]);
        /*
        foreach (var iatm in Enumerable.Range(0, Natm))
        {
            int p0 = mol.OffsetsByAtoms[iatm];
            int p1 = mol.OffsetsByAtoms[iatm + 1];
            NRange sa = p0..p1;
            NRange atmDim = new NRange(iatm * 3, iatm * 3 + 3);
            overlapDerivative[atmDim, sa, ..] = Tensor.Negate<double>(int1e_ipovlp[.., sa, ..]);
            hamiltonianDerivative[atmDim, sa, ..] = Tensor.Negate<double>(Tensor.Add<double>(int1e_ipkin, int1e_ipnucl))[.., sa, ..];

            using (mol.RinvAt(iatm))
            {
                var int1e_iprinv = mol.InvokeIntor("int1e_iprinv");
                hamiltonianDerivative[atmDim, .., ..] =
                    Tensor.Subtract<double>(
                        hamiltonianDerivative[atmDim, .., ..],
                        Tensor.Multiply<double>(int1e_iprinv, z_a[iatm]));
            }
        }
        Console.WriteLine(overlapDerivative.ToString());
        Console.WriteLine(hamiltonianDerivative.ToString());
        PrintRMS(overlapDerivative, nameof(overlapDerivative));
        PrintRMS(hamiltonianDerivative, nameof(hamiltonianDerivative));
        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < nao; j++)
            {
                for (int k = 0; k < nao; k++)
                {
                    overlapDerivative[i, j, k] += overlapDerivative[i, k, j];
                    hamiltonianDerivative[i, j, k] += hamiltonianDerivative[i, k, j];
                }
            }
        }
        */
        ///*
        for (var iatm = 0; iatm < Natm; iatm++)
        {
            var sa = mol.EnumerateByAtom(iatm);
            using (mol.RinvAt(iatm))
            {
                var int1e_iprinv = mol.InvokeIntor("int1e_iprinv");
                for (int idim = 0; idim < 3; idim++)
                {
                    int totalDim = iatm * 3 + idim;
                    for (var k = 0; k < nao; k++)
                    {
                        foreach (var j in sa)
                        {
                            double value = int1e_ipovlp[idim, j, k];
                            overlapDerivative[totalDim, j, k] -= value;
                            overlapDerivative[totalDim, k, j] -= value;
                            value = int1e_core[idim, j, k];
                            hamiltonianDerivative[totalDim, j, k] -= value;
                            hamiltonianDerivative[totalDim, k, j] -= value;
                        }
                        for (int j = 0; j < nao; j++)
                        {
                            double value = z_a[iatm] * int1e_iprinv[idim, j, k];
                            hamiltonianDerivative[totalDim, j, k] -= value;
                            hamiltonianDerivative[totalDim, k, j] -= value;
                        }
                    }
                }
            }
        }
        #endregion

        #region eriDerivative
        var int2e_ip1 = mol.InvokeIntor("int2e_ip1");
        eriDerivative = NDArray.Create([length, Nao, Nao, Nao, Nao]);
        for (var iatm = 0; iatm < Natm; iatm++)
        {
            var sa = mol.EnumerateByAtom(iatm);
            for (int idim = 0; idim < 3; idim++)
            {
                int totalDim = iatm * 3 + idim;
                for (var i = 0; i < nao; i++)
                {
                    for (var j = 0; j < nao; j++)
                    {
                        for (var k = 0; k < nao; k++)
                        {
                            foreach (var s in sa)
                            {
                                eriDerivative[totalDim, s, i, j, k] -= int2e_ip1[idim, s, i, j, k];
                                eriDerivative[totalDim, i, s, j, k] -= int2e_ip1[idim, s, i, j, k];
                                eriDerivative[totalDim, i, j, s, k] -= int2e_ip1[idim, s, k, i, j];
                                eriDerivative[totalDim, i, j, k, s] -= int2e_ip1[idim, s, k, i, j];
                            }
                        }
                    }
                }
            }
        }
        #endregion

        #region electronicEnergyDerivative
        overlapDerivativeMO = NDArray.Create([length, Nao, Nao]);
        var fock_mo = AO2MO(coeff, fock);
        Span<NRange> ranges = [NRange.All, NRange.All, NRange.All];
        overlapDerivativeMO[.., .., ..] = AO2MO(coeff, overlapDerivative);
        elecEnergyDerivative = NDArray.Create([length]);
        var sub1 = NDArray.Create([length]);
        var sub2 = NDArray.Create([length]);
        var sub3 = NDArray.Create([length]);
        var sub4 = NDArray.Create([length]);
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
        elecEnergyDerivative = elecEnergyDerivative + sub1;
        //sub2 = torch.einsum("duvkl,uv,kl -> d", eriDerivative, density, density);

        var uv = density;
        var kl = density;
        var dkl = NDArray.Create([Nao, Nao]);
        for (int d = 0; d < length; d++)
        {
            var duvkl = eriDerivative;
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
                    dkl[k, l] = idkl;
                }
            }
            double id = 0;
            for (int k = 0; k < Nao; k++)
            {
                for (int l = 0; l < Nao; l++)
                {
                    id += kl[k, l] * dkl[k, l];
                }
            }
            sub2[d] += id / 2;
        }
        elecEnergyDerivative = elecEnergyDerivative + sub2;
        //sub3 = -torch.einsum("dukvl,uv,kl -> d", eriDerivative, density, density);
        var dukvl = eriDerivative;
        dkl = NDArray.Create([length, Nao, Nao]);
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
            sub3[d] += id / 4;
        }
        elecEnergyDerivative = elecEnergyDerivative - sub3;
        //sub4 -= torch.einsum("dij, ij -> d", overlapDerivativeMO, fock_mo);
        //sub4 = np.einsum('atij, ij -> at', self.overlap_mo_derivative[:, :, so, so], fock_mo[so, so]) 
        var dij = overlapDerivativeMO;
        var ij = fock_mo;
        for (int d = 0; d < length; d++)
        {
            double id = 0;
            for (int i = 0; i < nocc; i++)
            {
                for (int j = 0; j < nocc; j++)
                {
                    id += dij[d, i, j] * ij[i, j];
                }
            }
            sub4[d] += id * 2;
        }
        elecEnergyDerivative = elecEnergyDerivative - sub4;
        PrintForce(elecEnergyDerivative);
        //elecEnergyDerivative = Vector2Tensor(sub1 + sub2 / 2 - sub3 / 4 - sub4 * 2);
        #endregion

        #region nuclear repulsion derivative
        nuclearRepulsionDerivative = NDArray.Create([length]);
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
                    var coorSpan = coor.AsSpan();
                    var z_ij = z_i * z_j;
                    value -= z_ij * coorSpan[idim] / Math.Pow(coor.Length, 3);
                }

                nuclearRepulsionDerivative[totalDim] = value;
            }
        }
        #endregion

        for (int i = 0; i < Natm * 3; i++)
        {
            force[i] = elecEnergyDerivative[i] + nuclearRepulsionDerivative[i];
        }
        PrintForce(force);
    }

    private static NDArray AO2MO(NDArray moCoeffs, NDArray getter)
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
        for (var idim = 0; idim < dims.Length; idim++)
        {
            var dim = dims[idim];
            if (dims.Length == 3 && idim == 0)
            {
                continue;
            }
            if (dim != length)
            {
                throw new ArgumentException("The input matrix must have the same size as the moCoeffs.");
            }
        }
        NDArray temp1 = NDArray.Create(dims);
        NDArray temp2 = NDArray.Create(dims);

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
            case 3:
                {

                    var duv = getter;
                    var up = moCoeffs;
                    var dpv = temp1;
                    for (nint idim = 0; idim < getter.Lengths[0]; idim++)
                    {
                        for (int p = 0; p < length; p++)
                        {
                            for (int v = 0; v < length; v++)
                            {
                                var idpv = 0.0;
                                for (int u = 0; u < length; u++)
                                {
                                    idpv += duv[idim, u, v] * up[u, p];
                                }
                                dpv[idim, p, v] = idpv;
                            }
                        }
                    }
                    var vq = moCoeffs;
                    var dpq = temp2;
                    for (nint idim = 0; idim < getter.Lengths[0]; idim++)
                    {
                        for (int p = 0; p < length; p++)
                        {
                            for (int q = 0; q < length; q++)
                            {
                                var idpq = 0.0;
                                for (int v = 0; v < length; v++)
                                {
                                    idpq += dpv[idim, p, v] * vq[v, q];
                                }
                                dpq[idim, p, q] = idpq;
                            }
                        }
                    }
                    return dpq;
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

    private void PrintForce(NDArray force)
    {
        StringBuilder sb = new();
        for (int i = 0; i < Natm; i++)
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
