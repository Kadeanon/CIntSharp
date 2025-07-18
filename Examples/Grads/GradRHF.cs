using CintSharp.DataStructures;
using CintSharp.Intor;
using Examples.Methods;
using SimpleHelpers;
using SimpleHelpers.MultiAlg;
using SimpleHelpers.MultiAlg.TensorContract;
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

namespace Examples.Grads;

public class GradRHF
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


    public NDArray? OverlapDerivative { get; private set; }
    public NDArray? HamiltonianDerivative { get; private set; }
    public NDArray? EriDerivative { get; private set; }
    public NDArray? ElecEnergyDerivative { get; private set; }
    public NDArray? NuclearRepulsionDerivative { get; private set; }
    public NDArray? OverlapDerivativeMO { get; private set; }

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

        var int1e_ipcore = int1e_ipkin.AddedBy(int1e_ipnucl);


        #endregion

        #region overlapDerivativ and hamiltonianDerivative
        //derivatives in AO basis
        OverlapDerivative = NDArray.Create([length, nao, nao]);
        HamiltonianDerivative = NDArray.Create([length, nao, nao]);
        for (var iatm = 0; iatm < Natm; iatm++)
        {
            var sa = mol.EnumerateByAtom(iatm);
            var sr = mol.RangesByAtoms[iatm];
            using var _ = mol.RinvAt(iatm);
            var int1e_iprinv = mol.InvokeIntor("int1e_iprinv");
            for (int idim = 0; idim < 3; idim++)
            {
                int totalDim = iatm * 3 + idim;
                for (var k = 0; k < nao; k++)
                {
                    foreach (var j in sa)
                    {
                        double value = int1e_ipovlp[idim, j, k];
                        OverlapDerivative[totalDim, j, k] -= value;
                        OverlapDerivative[totalDim, k, j] -= value;
                        value = int1e_ipcore[idim, j, k];
                        HamiltonianDerivative[totalDim, j, k] -= value;
                        HamiltonianDerivative[totalDim, k, j] -= value;
                    }
                    for (int j = 0; j < nao; j++)
                    {
                        double value = z_a[iatm] * int1e_iprinv[idim, j, k];
                        HamiltonianDerivative[totalDim, j, k] -= value;
                        HamiltonianDerivative[totalDim, k, j] -= value;
                    }
                }
            }
        }
        #endregion

        #region eriDerivative
        var int2e_ip1 = mol.InvokeIntor("int2e_ip1");
        EriDerivative = NDArray.Create([Natm, 3, Nao, Nao, Nao, Nao]);
        for (var iatm = 0; iatm < Natm; iatm++)
        {
            var sr = mol.RangesByAtoms[iatm];
            var int2e_ip1Slice = int2e_ip1[.., sr];
            EriDerivative[iatm, .., sr, .., .., ..]
                -= int2e_ip1Slice;
        }
        for (var iatm = 0; iatm < Natm; iatm++)
            {
                var sr = mol.RangesByAtoms[iatm];
                var int2e_ip1Slice = int2e_ip1[.., sr];
                EriDerivative[iatm, .., .., sr, .., ..]
                -= int2e_ip1Slice.Transpose([0, 2, 1, 3, 4]);
        }
        for (var iatm = 0; iatm < Natm; iatm++)
        {
            var sr = mol.RangesByAtoms[iatm];
            var int2e_ip1Slice = int2e_ip1[.., sr];
            EriDerivative[iatm, .., .., .., sr, ..]
                -= int2e_ip1Slice.Transpose([0, 3, 4, 1, 2]);
        }
        for (var iatm = 0; iatm < Natm; iatm++)
        {
            var sr = mol.RangesByAtoms[iatm];
            var so2 = mol.EnumerateByAtom(iatm);
            foreach(var idim in Enumerable.Range(0, 3))
                foreach (var i in Enumerable.Range(0, nao))
                    foreach (var j in Enumerable.Range(0, nao))
                        foreach (var k in Enumerable.Range(0, nao))
                            foreach (var s in so2)
                                EriDerivative[iatm, idim, i, j, k, s] -=
                                    int2e_ip1[idim, s, k, i, j];
        }
        var eriDerivativeOld = EriDerivative;
        EriDerivative = EriDerivative.Reshape([-1, Nao, Nao, Nao, Nao]);
        #endregion

        #region electronicEnergyDerivative
        OverlapDerivativeMO = NDArray.Create([length, Nao, Nao]);
        var fock_mo = AO2MO(fock, coeff[.., ..nocc]);
        OverlapDerivativeMO[.., .., ..] = AO2MO(OverlapDerivative, coeff);
        ElecEnergyDerivative = NDArray.Contract("tij,ij->t",
            HamiltonianDerivative, density);
        NDArray.Einsum(ElecEnergyDerivative, "kl,duvkl,uv->d",
            0.5, density, EriDerivative, density);
        NDArray.Einsum(ElecEnergyDerivative, "kl,dukvl,uv->d",
            -0.25, density, EriDerivative, density);
        NDArray.Contract(ElecEnergyDerivative, "dij,ij->d",
            -2, OverlapDerivativeMO[.., ..nocc, ..nocc], 
            fock_mo[..nocc, ..nocc]);
        #endregion

        #region nuclear repulsion derivative
        NuclearRepulsionDerivative = NDArray.Create([length]);
        for (var totalDim = 0; totalDim < Natm * 3; totalDim++)
        {
            var iatm = totalDim / 3;
            var z_i = z_a[iatm];
            var coor_i = ac[iatm];
            double value = 0;
            for (var jatm = 0; jatm < Natm; jatm++)
            {
                if (jatm == iatm)
                    continue;

                var z_j = z_a[jatm];
                var coor_j = ac[jatm];
                var coor = coor_i - coor_j;
                var z_ij = z_i * z_j;
                value -= z_ij * coor[totalDim % 3] / Math.Pow(coor.Length, 3);
            }
            NuclearRepulsionDerivative[totalDim] = value;
        }
        #endregion

        force = ElecEnergyDerivative + NuclearRepulsionDerivative;
    }

    private static NDArray AO2MO(NDArray ao, NDArray coeffs)
    {
        var dims = ao.Lengths;
        if (coeffs.Rank != 2)
            throw new ArgumentException("The input matrix must be 2D.");
        nint length = coeffs.Lengths[0];
        for (var idim = 0; idim < dims.Length; idim++)
        {
            var dim = dims[idim];
            if (dims.Length == 3 && idim == 0)
                continue;
            if (dim != length)
                throw new ArgumentException(
                    "The input matrix must have the same size as the moCoeffs.");
        }

        return dims.Length switch
        {
            2 => NDArray.Einsum("uv,up,vq->pq",
                                    ao, coeffs, coeffs),
            3 => NDArray.Einsum("duv,up,vq->dpq",
                                    ao, coeffs, coeffs),
            4 => NDArray.Einsum("uvkl,up,vq,kr,ls->pqrs",
                                    ao, coeffs, coeffs, coeffs, coeffs),
            _ => throw new ArgumentException("The input matrix must be 2D or 4D."),
        };
    }
}
