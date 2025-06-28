using SimpleHelpers;
using SimpleHelpers.MultiAlg;
using static SimpleHelpers.MultiAlg.TensorContract.ContractMethods;
using System.Buffers;

namespace Examples
{
    internal class CCD
    {
        public nint nocc;
        public nint nvir;
        public nint nmo;
        public NDArray mo_energy;
        public NDArray eris;
        public int max_cycle;
        public double tol;
        public NDArray t2;
        public double eccd;
        NRange occRange;
        NRange virRange;

        public bool Run(int max_cycle = 50, double tol = 1e-8)
        {
            var eris_oovv = eris[occRange, occRange, virRange, virRange];
            var e_init = Energy(eris_oovv);
            Console.WriteLine($"CCD initial correlation energy: {e_init:F12}");
            bool converged = false;
            eccd = 0.0;
            int i = 0;
            for (; i < max_cycle; i++)
            {
                var t2new = UpdateAmps(t2, eris, mo_energy);
                var norm_residual = (t2new - t2).Norm();
                t2 = t2new;
                eccd = Energy(eris_oovv);
                Console.WriteLine($"CCD iter {i}: " +
                    $"E_corr = {eccd:F12}, residual = {norm_residual:F12}");
                if (norm_residual < tol)
                {
                    converged = true;
                    break;
                }
            }
            if (!converged)
                Console.WriteLine("CCD did not converge within the maximum number of iterations.");
            else
                Console.WriteLine($"CCD converged in {i + 1} iterations with correlation energy: {eccd:F12}");

            return converged;
        }

        public NDArray InitAmps(NDArray eris_vvoo, NDArray mo_ene)
        {
            NDArray t2 = NDArray.CreateUninitialized([nvir, nvir, nocc, nocc]);
            for (int a = 0; a < nvir; a++)
            {
                var mo_a = mo_ene[nocc + a];
                for (int b = 0; b < nvir; b++)
                {
                    var mo_b = mo_ene[nocc + b];
                    var mo_ab = mo_a + mo_b;
                    for (int i = 0; i < nocc; i++)
                    {
                        var mo_i = mo_ene[i];
                        for (int j = 0; j < nocc; j++)
                        {
                            t2[a, b, i, j] =
                                eris_vvoo[a, b, i, j] / (mo_i + mo_ene[j] - mo_ab);
                        }
                    }
                }
            }
            return t2;
        }

        public NDArray UpdateAmps(NDArray t2, NDArray eris, NDArray mo_ene)
        {
            var eris_oooo = eris[occRange, occRange, occRange, occRange];
            var eris_oovv = eris[occRange, occRange, virRange, virRange];
            var eris_vvoo = eris[virRange, virRange, occRange, occRange];
            var eris_vvvv = eris[virRange, virRange, virRange, virRange];
            var eris_ovvo = eris[occRange, virRange, virRange, occRange];

            var ft_vv = Contract(
                "afmn,mnbf->ab", -0.5, t2, eris_oovv);
            var ft_oo = Contract(
                "efin,jnef->ij", 0.5, t2, eris_oovv);
            var w_oooo = eris_oooo + Contract(
                "efij,mnef->mnij", 0.25, t2, eris_oovv);
            var w_vvvv = eris_vvvv + Contract(
                "abmn,mnef->abef", 0.25, t2, eris_oovv);
            var w_ovvo = eris_ovvo - Contract(
                "fbjn,mnef->mbej", 0.5, t2, eris_oovv);
            var tmp = Contract(
                "aeij,be->abij", t2, ft_vv);
            var t2new = eris_vvoo + tmp - tmp.Transpose(1, 0, 2, 3);
            tmp = Contract(
                "abim,mj->abij", t2, ft_oo);
            t2new
                .SubtractedBy(tmp)
                .AddedBy(tmp.Transpose(0, 1, 3, 2))
                .AddedBy(Contract(
                    "abmn,mnij->abij", 0.5, t2, w_oooo))
                .AddedBy(Contract(
                    "efij,abef->abij", 0.5, t2, w_vvvv));
            tmp = Contract(
                "aeim,mbej->abij", t2, w_ovvo);
            t2new
                .AddedBy(tmp)
                .SubtractedBy(tmp.Transpose(1, 0, 2, 3))
                .SubtractedBy(tmp.Transpose(0, 1, 3, 2))
                .AddedBy(tmp.Transpose(1, 0, 3, 2));
            var mo_occ = mo_ene[..(int)nocc];
            var mo_vir = mo_ene[(int)nocc..];
            for (int i = 0; i < nocc; i++)
                for (int j = 0; j < nocc; j++)
                    for (int a = 0; a < nvir; a++)
                        for (int b = 0; b < nvir; b++)
                            t2new[a, b, i, j] /= mo_occ[i] + mo_occ[j] - mo_vir[a] - mo_vir[b];

            return t2new;
        }

        public double Energy(NDArray eris_oovv)
        {
            double e = 0.0;
            for (int i = 0; i < nocc; i++)
            {
                for (int j = 0; j < nocc; j++)
                {
                    for (int a = 0; a < nvir; a++)
                    {
                        for (int b = 0; b < nvir; b++)
                        {
                            e += t2[a, b, i, j] * eris_oovv[i, j, a, b];
                        }
                    }
                }
            }
            return e / 4;
        }

        public static NDArray AO2MO(NDArray ao, NDArray coeff)
        {
            if (ao.Rank != 4)
                throw new ArgumentException("AO integral must be a 4D array.");
            var ajkl = Contract
                ("abcd,ai->ibcd", ao, coeff);
            var ijcd = Contract
                ("ibcd,bj->ijcd", ajkl, coeff);
            var ijkd = Contract
                ("ijcd,ck->ijkd", ijcd, coeff);
            var ijkl = Contract
                ("ijkd,dl->ijkl", ijkd, coeff);
            return ijkl;
        }

        public static NDArray AO2MO(NDArray ao, NDArray ai, NDArray bj, NDArray ck, NDArray dl)
        {
            if (ao.Rank != 4)
                throw new ArgumentException("AO integral must be a 4D array.");
            var ajkl = SimpleContract
                ("abcd,ai->ibcd", ao, ai);
            var ijcd = SimpleContract
                ("ibcd,bj->ijcd", ajkl, bj);
            var ijkd = SimpleContract
                ("ijcd,ck->ijkd", ijcd, ck);
            var ijkl = SimpleContract
                ("ijkd,dl->ijkl", ijkd, dl);
            return ijkl;
        }

        public NDArray MakeErisIncore(RHF mf)
        {
            var mo_coeff = mf.C!.AsNDArray();
            var nao = mo_coeff.Lengths[0];
            var ao_eri = mf.ERI;
            var mo_a = NDArray.Create([nao, nmo]);
            mo_a.Reshape([nao, nao, 2])[.., .., 0] = mo_coeff;
            var mo_b = NDArray.Create([nao, nmo]);
            mo_b.Reshape([nao, nao, 2])[.., .., 1] = mo_coeff;
            var eri = AO2MO(ao_eri, mo_a);
            eri += AO2MO(ao_eri, mo_b);
            var eri1 = AO2MO(ao_eri, mo_a, mo_a, mo_b, mo_b);
            eri += eri1;
            eri += eri1.Transpose();
            return eri;
        }

        public CCD(RHF mf)
        {
            var nao = mf.Envs.NAO;
            mo_energy = NDArray.CreateUninitialized([nao * 2]);
            var unfold = mo_energy.Reshape([nao, 2]);
            var ao_ene = mf.Es!.AsNDArray();
            unfold[.., 1] = ao_ene;
            unfold[.., 0] = ao_ene;
            nocc = mf.Atoms.Sum(atom => atom.AtomNumber);
            nmo = nao * 2;
            nvir = nmo - nocc;
            Console.WriteLine($"Nocc={nocc}, nvir={nvir}");
            eris = MakeErisIncore(mf);
            eris = eris.Transpose(0, 2, 1, 3) - eris.Transpose(0, 2, 3, 1);
            max_cycle = 50;
            tol = 1e-8;
            occRange = new NRange(0, nocc);
            virRange = new NRange(nocc, nmo);
            var eris_vvoo = eris[virRange, virRange, occRange, occRange];
            t2 = InitAmps(eris_vvoo, mo_energy);
        }
    }
}
