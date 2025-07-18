using SimpleHelpers;
using SimpleHelpers.MultiAlg;
using ND = SimpleHelpers.MultiAlg.NDArray;
using System.Buffers;
using SimpleHelpers.LinearAlg;

namespace Examples.Methods
{
    public class CCD
    {
        public int MaxCycles { get; set; }
        public double Tol { get; set; }
        public nint NOcc { get; }
        public nint NVir { get; }
        public nint NMO { get; }
        public Vector MOEnergies { get; }
        public NDArray Eris { get; }
        public NDArray T2 { get; private set; }
        public double Eccd { get; private set; }
        NRange OccRange { get; }
        NRange VirRange { get; }

        public CCD(RHF mf)
        {
            var nao = mf.Envs.NAO;
            MOEnergies = Vector.Create(nao * 2, uninited: true);
            var unfold = MOEnergies.AsNDArray().Reshape([nao, 2]);
            var aoEnergies = mf.Es!.AsNDArray();
            unfold[.., 1] = aoEnergies;
            unfold[.., 0] = aoEnergies;
            NOcc = mf.Atoms.Sum(atom => atom.AtomNumber);
            NMO = nao * 2;
            NVir = NMO - NOcc;
            Console.WriteLine($"Nocc={NOcc}, nvir={NVir}");
            Eris = ToMOERI(mf);
            Eris = Eris.Transpose(0, 2, 1, 3) - Eris.Transpose(0, 2, 3, 1);
            MaxCycles = 50;
            Tol = 1e-8;
            OccRange = new NRange(0, NOcc);
            VirRange = new NRange(NOcc, NMO);
            var eris_vvoo = Eris[VirRange, VirRange, OccRange, OccRange];
            T2 = InitAmps(eris_vvoo, MOEnergies);
            Eccd = 0.0;
        }

        public bool Run(int max_cycle = 50, double tol = 1e-8)
        {
            var eris_oovv = Eris[OccRange, OccRange, VirRange, VirRange];
            var e_init = Energy(eris_oovv);
            Console.WriteLine($"CCD initial correlation energy: {e_init:F12}");
            bool converged = false;
            int i = 0;
            for (; i < max_cycle; i++)
            {
                var t2new = UpdateAmps(T2, Eris, MOEnergies);
                var norm_residual = (t2new - T2).Norm();
                T2 = t2new;
                Eccd = Energy(eris_oovv);
                Console.WriteLine($"CCD iter {i}: " +
                    $"E_corr = {Eccd:F12}, residual = {norm_residual:F12}");
                if (norm_residual < tol)
                {
                    converged = true;
                    break;
                }
            }
            if (!converged)
                Console.WriteLine("CCD did not converge within the maximum number of iterations.");
            else
                Console.WriteLine($"CCD converged in {i + 1} iterations with correlation energy: {Eccd:F12}");

            return converged;
        }

        public NDArray InitAmps(NDArray eris_vvoo, Vector mo_ene)
        {
            NDArray t2 = NDArray.CreateUninitialized([NVir, NVir, NOcc, NOcc]);
            var mo_occ = mo_ene[OccRange];
            var mo_vir = mo_ene[VirRange];
            for (int a = 0; a < NVir; a++)
            {
                for (int b = 0; b < NVir; b++)
                {
                    var mo_ab = mo_vir[a] + mo_vir[b];
                    for (int i = 0; i < NOcc; i++)
                    {
                        var mo_i = mo_occ[i];
                        for (int j = 0; j < NOcc; j++)
                        {
                            var mo_ij = mo_i + mo_occ[i];
                            t2[a, b, i, j] =
                                eris_vvoo[a, b, i, j] / (mo_ij - mo_ab);
                        }
                    }
                }
            }
            return t2;
        }

        public NDArray UpdateAmps(NDArray t2, NDArray eris, Vector mo_ene)
        {
            var residual = ComputeResidual(t2, eris);
            var mo_occ = mo_ene[OccRange];
            var mo_vir = mo_ene[VirRange];
            for (int a = 0; a < NVir; a++)
            {
                for (int b = 0; b < NVir; b++)
                {
                    var mo_ab = mo_vir[a] + mo_vir[b];
                    for (int i = 0; i < NOcc; i++)
                    {
                        var mo_i = mo_occ[i];
                        for (int j = 0; j < NOcc; j++)
                        {
                            var mo_ij = mo_i + mo_occ[j];
                            var d_abij = mo_ij - mo_ab;
                            residual[a, b, i, j] /= d_abij;
                        }
                    }
                }
            }

            return residual;
        }

        public NDArray ComputeResidual(NDArray t2, NDArray eris)
        {
            var eris_oooo = eris[OccRange, OccRange, OccRange, OccRange];
            var eris_oovv = eris[OccRange, OccRange, VirRange, VirRange];
            var eris_vvoo = eris[VirRange, VirRange, OccRange, OccRange];
            var eris_vvvv = eris[VirRange, VirRange, VirRange, VirRange];
            var eris_voov = eris[VirRange, OccRange, OccRange, VirRange];

            var residual = eris_vvoo.Clone();
            ND.Contract(residual, "abcd,cdij->abij", 0.5, eris_vvvv, t2);
            ND.Contract(residual, "klij,abkl->abij", 0.5, eris_oooo, t2);
            Pabij(residual, "akic,bcjk->abij", +1.0, eris_voov, t2);
            Pab(residual, "klcd,bdkl,acij->abij", -0.5, eris_oovv, t2, t2);
            Pij(residual, "klcd,cdjl,abik->abij", -0.5, eris_oovv, t2, t2);
            ND.Einsum(residual, "klcd,cdij,abkl->abij", 0.25, eris_oovv, t2, t2);
            Pij(residual, "klcd,acik,bdjl->abij", 1.0, eris_oovv, t2, t2);
            return residual;
        }

        private static void Pab(NDArray self, string expression,
            double scalar, params ReadOnlySpan<NDArray> tensors)
        {
            var result = ND.Einsum(expression, scalar, tensors);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([1, 0, 2, 3]));
        }

        private static void Pij(NDArray self, string expression,
            double scalar, params ReadOnlySpan<NDArray> tensors)
        {
            var result = ND.Einsum(expression, scalar, tensors);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([0, 1, 3, 2]));
        }

        private static void Pabij(NDArray self, string expression,
            double scalar, NDArray left, NDArray right)
        {
            var result = ND.Contract(expression, scalar, left, right);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([1, 0, 2, 3]))
                .SubtractedBy(result.Transpose([0, 1, 3, 2]))
                .AddedBy(result.Transpose([1, 0, 3, 2]));
        }

        public double Energy(NDArray eris_oovv)
        {
            double e = 0.0;
            for (int i = 0; i < NOcc; i++)
            {
                for (int j = 0; j < NOcc; j++)
                {
                    for (int a = 0; a < NVir; a++)
                    {
                        for (int b = 0; b < NVir; b++)
                        {
                            e += T2[a, b, i, j] * eris_oovv[i, j, a, b];
                        }
                    }
                }
            }
            return e / 4;
        }

        public static NDArray AO2MO(NDArray ao, NDArray coeff)
            => ND.Einsum("pqrs,pi,qj,rk,sl->ijkl",
                ao, coeff, coeff, coeff, coeff);

        public static NDArray AO2MO(NDArray ao,
            NDArray pi, NDArray qj, NDArray rk, NDArray sl)
            => ND.Einsum("pqrs,pi,qj,rk,sl->ijkl",
                ao, pi, qj, rk, sl);

        public NDArray ToMOERI(RHF mf)
        {
            var mo_coeff = mf.C!.AsNDArray();
            var nao = mo_coeff.Lengths[0];
            var ao_eri = mf.ERI;
            var mo_a = NDArray.Create([nao, NMO]);
            mo_a.Reshape([nao, nao, 2])[.., .., 0] = mo_coeff;
            var mo_b = NDArray.Create([nao, NMO]);
            mo_b.Reshape([nao, nao, 2])[.., .., 1] = mo_coeff;
            var eri = AO2MO(ao_eri, mo_a);
            eri += AO2MO(ao_eri, mo_b);
            var eri1 = AO2MO(ao_eri, mo_a, mo_a, mo_b, mo_b);
            eri += eri1;
            eri += eri1.Transpose();
            return eri;
        }
    }
}
