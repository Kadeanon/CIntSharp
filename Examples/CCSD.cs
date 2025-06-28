using SimpleHelpers;
using SimpleHelpers.MultiAlg;
using static SimpleHelpers.MultiAlg.TensorContract.ContractMethods;
using System.Buffers;
using System.Diagnostics.CodeAnalysis;

namespace Examples
{
    internal class CCSD
    {
        public nint nocc;
        public nint nvir;
        public nint nmo;
        public NDArray mo_energy;
        public NDArray eris;
        public int max_cycle;
        public double tol;
        public NDArray t1;
        public NDArray t2;
        public double eccd;
        NRange occRange;
        NRange virRange;

        public bool Run(int max_cycle = 50, double tol = 1e-8)
        {
            var mo_ene = mo_energy;
            var eris_oovv = eris[occRange, occRange, virRange, virRange];
            var e_init = Energy(t1, t2, eris_oovv);
            Console.WriteLine($"CCSD initial correlation energy: {e_init:F12}");
            bool converged = false;
            eccd = 0.0;
            int i = 0;
            for (; i < max_cycle; i++)
            {
                var (t1new, t2new) = UpdateAmps(eris, mo_ene);
                var norm_residual = (t2new - t2).Norm();
                t1 = t1new;
                t2 = t2new;
                eccd = Energy(t1new, t2new, eris_oovv);
                Console.WriteLine($"CCSD iter {i}: " +
                    $"E_corr = {eccd:F12}, residual = {norm_residual:F12}");
                if (norm_residual < tol)
                {
                    converged = true;
                    break;
                }
            }
            if (!converged)
                Console.WriteLine("CCSD did not converge within the maximum number of iterations.");
            else
                Console.WriteLine($"CCSD converged in {i + 1} iterations with correlation energy: {eccd:F12}");

            return converged;
        }

        [MemberNotNull(nameof(t1), nameof(t2))]
        public void InitAmps(NDArray eris_vvoo, NDArray mo_ene)
        {
            t2 = NDArray.CreateUninitialized([nvir, nvir, nocc, nocc]);
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
            t1 = NDArray.Create([nvir, nocc]);
        }

        public(NDArray t1, NDArray t2)
            UpdateAmps(NDArray eris, NDArray mo_ene)
        {
            var eris_oooo = eris[occRange, occRange, occRange, occRange];
            var eris_oovv = eris[occRange, occRange, virRange, virRange];
            var eris_vvoo = eris[virRange, virRange, occRange, occRange];
            var eris_vvvv = eris[virRange, virRange, virRange, virRange];
            var eris_ovvo = eris[occRange, virRange, virRange, occRange];
            var eris_ovvv = eris[occRange, virRange, virRange, virRange];
            var eris_ooov = eris[occRange, occRange, occRange, virRange];
            var eris_vvvo = eris[virRange, virRange, virRange, occRange];
            var eris_ovoo = eris[occRange, virRange, occRange, occRange];

            var tmp = Contract("ai,bj->abij", t1, t1);
            tmp -= tmp.Transpose(1, 0, 2, 3);
            var tau = t2 + tmp;
            var taut = t2 + 0.5 * tmp;
            var ft_vv = Contract("fm,mafb->ab", t1, eris_ovvv) -
                0.5 * Contract("afmn,mnbf->ab", taut, eris_oovv);
            var ft_oo = Contract("en,inje->ij", t1, eris_ooov) +
                0.5 * Contract("efjn,inef->ij", taut, eris_oovv);
            var ft_ov = Contract("fn,inaf->ia", t1, eris_oovv);
            tmp = Contract("ej,mnie->mnij", t1, eris_ooov);
            tmp -= tmp.Transpose(0, 1, 3, 2);
            var w_oooo = eris_oooo +
                0.25 * Contract("efij,mnef->mnij", tau, eris_oovv) + tmp;
            tmp = Contract("bm,maef->abef", t1, eris_ovvv);
            tmp -= tmp.Transpose(1, 0, 2, 3);
            var w_vvvv = eris_vvvv +
                0.25 * Contract("abmn,mnef->abef", tau, eris_oovv) + tmp;
            var w_ovvo = eris_ovvo +
                Contract("fj,mbef->mbej", t1, eris_ovvv) +
                Contract("bn,mnje->mbej", t1, eris_ooov);
            tmp = 0.5 * t2 + Contract("fj,bn->fbjn", t1, t1);
            w_ovvo -= Contract("fbjn,mnef->mbej", tmp, eris_oovv);

            var t1new = NDArray.ZeroLike(t1);
            t1new += Contract("ei,ae->ai", t1, ft_vv);
            t1new -= Contract("am,mi->ai", t1, ft_oo);
            t1new += Contract("aeim,me->ai", t2, ft_ov);
            t1new += Contract("fn,nafi->ai", t1, eris_ovvo);
            t1new -= 0.5 * Contract("efim,maef->ai", t2, eris_ovvv);
            t1new -= 0.5 * Contract("aemn,mnie->ai", t2, eris_ooov);
            for (int i = 0; i < nocc; i++)
            {
                for (int a = 0; a < nvir; a++)
                {
                    t1new[a, i] /= (mo_ene[i] - mo_ene[nocc + a]);
                }
            }


            ft_vv -= 0.5 * Contract("bm,me->be", t1, ft_ov);
            tmp = Contract("aeij,be->abij", t2, ft_vv);
            var t2new = eris_vvoo + tmp - tmp.Transpose(1, 0, 2, 3);
            ft_oo += 0.5 * Contract("ej,me->mj", t1, ft_ov);
            tmp = Contract("abim,mj->abij", t2, ft_oo);
            t2new -= tmp - tmp.Transpose(0, 1, 3, 2);
            t2new += 0.5 * Contract("abmn,mnij->abij", tau, w_oooo);
            t2new += 0.5 * Contract("efij,abef->abij", tau, w_vvvv);
            tmp = Contract("am,mbej->abej", t1, eris_ovvo);
            tmp = Contract("ei,abej->abij", t1, tmp);
            tmp = Contract("aeim,mbej->abij", t2, w_ovvo) - tmp;
            t2new += tmp - tmp.Transpose(1, 0, 2, 3)
                - tmp.Transpose(0, 1, 3, 2) + tmp.Transpose(1, 0, 3, 2);
            tmp = Contract("ei,abej->abij", t1, eris_vvvo);
            t2new += tmp - tmp.Transpose(0, 1, 3, 2);
            tmp = Contract("am,mbij->abij", t1, eris_ovoo);
            t2new -= tmp - tmp.Transpose(1, 0, 2, 3);

            var mo_occ = mo_ene[..(int)nocc];
            var mo_vir = mo_ene[(int)nocc..];
            for (int i = 0; i < nocc; i++)
                for (int j = 0; j < nocc; j++)
                    for (int a = 0; a < nvir; a++)
                        for (int b = 0; b < nvir; b++)
                            t2new[a, b, i, j] /= mo_occ[i] + mo_occ[j] - mo_vir[a] - mo_vir[b];

            return (t1new, t2new);
        }

        public double Energy(NDArray t1, NDArray t2, NDArray eris_oovv)
        {
            double e_t2 = 0.0;
            double e_t1 = 0.0;
            for (int i = 0; i < nocc; i++)
            {
                for (int j = 0; j < nocc; j++)
                {
                    for (int a = 0; a < nvir; a++)
                    {
                        for (int b = 0; b < nvir; b++)
                        {
                            var ijab = eris_oovv[i, j, a, b];
                            e_t2 += t2[a, b, i, j] * ijab;
                            e_t1 += t1[a, i] * t1[b, j] * ijab;
                        }
                    }
                }
            }
            return e_t2 / 4 + e_t1 / 2;
        }

        public static NDArray AO2MO(NDArray ao, NDArray coeff)
        {
            if (ao.Rank != 4)
                throw new ArgumentException("AO integral must be a 4D array.");
            GC.Collect();
            var ajkl = SimpleContract
                ("abcd,ai->ibcd", ao, coeff);
            var ijcd = SimpleContract
                ("ibcd,bj->ijcd", ajkl, coeff);
            var ijkd = SimpleContract
                ("ijcd,ck->ijkd", ijcd, coeff);
            var ijkl = SimpleContract
                ("ijkd,dl->ijkl", ijkd, coeff);
            return ijkl;
        }

        public static NDArray AO2MO(NDArray ao, NDArray ai, NDArray bj, NDArray ck, NDArray dl)
        {
            if (ao.Rank != 4)
                throw new ArgumentException("AO integral must be a 4D array.");
            GC.Collect();
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

        public CCSD(CCD ccd)
        {
            nocc = ccd.nocc;
            nvir = ccd.nvir;
            Console.WriteLine($"Nocc={nocc}, nvir={nvir}");
            nmo = ccd.nmo;
            mo_energy = ccd.mo_energy;
            eris = ccd.eris;
            max_cycle = ccd.max_cycle;
            tol = ccd.tol;
            t1 = NDArray.Create([nvir, nocc]);
            t2 = ccd.t2;
            eccd = ccd.eccd;
            occRange = new NRange(0, nocc);
            virRange = new NRange(nocc, nmo);
        }

        public CCSD(RHF mf)
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
            var eris_oovv = eris[occRange, occRange, virRange, virRange];
            var eris_vvoo = eris[virRange, virRange, occRange, occRange];
            InitAmps(eris_vvoo, mo_energy);
        }
    }
}
