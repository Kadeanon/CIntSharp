using SimpleHelpers;
using SimpleHelpers.MultiAlg;
using System.Buffers;
using System.Diagnostics.CodeAnalysis;
using SimpleHelpers.LinearAlg;

namespace Examples.Methods
{
    public class CCSD
    {
        public int MaxCycle { get; set; }

        public double Tol { get; set; }

        public nint NOcc { get; }

        public nint NVir { get; }

        public nint NMO { get; }

        public Vector MOEnergy { get; }

        public NDArray Eris { get; }

        public NDArray T1 { get; private set; }

        public NDArray T2 { get; private set; }

        //public CCSDDIIS DIIS { get; set; }

        public double Eccsd { get; private set; }

        NRange OccRange { get; }

        NRange VirRange { get; }

        public double[]? tempArray;

        public double[]? tempArray2;

        public double[]? tempArray3;

        public bool Run()
        {
            var start = DateTime.Now;
            var moEnergy = MOEnergy;
            var eris_vvoo = Eris[VirRange, VirRange, OccRange, OccRange];
            var e0 = Energy(T1, T2, eris_vvoo);
            Console.WriteLine($"CCSD initial correlation energy: {e0:F12}");
            bool converged = false;
            Eccsd = 0.0;
            int i = 0;
            for (; i < MaxCycle; i++)
            {
                UpdateAmps(Eris, moEnergy, out var t1new, out var t2new);
                var norm_residual = (t2new - T2).Norm();
                T1 = t1new;
                T2 = t2new;
                Eccsd = Energy(T1, T2, eris_vvoo);
                Console.WriteLine($"CCSD iter {i}: " +
                    $"E_corr = {Eccsd:F12}, Norm_residual = {norm_residual:F12}");
                if (norm_residual < Tol)
                {
                    converged = true;
                    break;
                }
            }
            if (!converged)
                Console.WriteLine("CCSD did not converge within the maximum number of iterations.");
            else
                Console.WriteLine($"CCSD converged in {i + 1} iterations with correlation energy: {Eccsd:F12}");

            var end = DateTime.Now;
            var span = end - start;
            Console.WriteLine($"CCSD calculation took {span.TotalSeconds:F2} seconds.");
            return converged;
        }

        [MemberNotNull(nameof(T1), nameof(T2))]
        public void InitAmps(NDArray eris_vvoo, Vector mo_ene)
        {
            T2 = NDArray.CreateUninitialized([NVir, NVir, NOcc, NOcc]);
            for (int a = 0; a < NVir; a++)
            {
                var mo_a = mo_ene[NOcc + a];
                for (int b = 0; b < NVir; b++)
                {
                    var mo_b = mo_ene[NOcc + b];
                    var mo_ab = mo_a + mo_b;
                    for (int i = 0; i < NOcc; i++)
                    {
                        var mo_i = mo_ene[i];
                        for (int j = 0; j < NOcc; j++)
                        {
                            T2[a, b, i, j] =
                                eris_vvoo[a, b, i, j] / (mo_i + mo_ene[j] - mo_ab);
                        }
                    }
                }
            }
            T1 = NDArray.Create([NVir, NOcc]);
        }

        public void UpdateAmps(NDArray eris, Vector mo_ene,
            out NDArray t1, out NDArray t2)
        {
            ComputeResidualMin(eris,
                out var t1Residual, out var t2Residual);

            t1 = t1Residual;
            var mo_occ = mo_ene[OccRange];
            var mo_vir = mo_ene[VirRange];
            for (int a = 0; a < NVir; a++)
                for (int i = 0; i < NOcc; i++)
                    t1[a, i] /= mo_occ[i] - mo_vir[a];

            t2 = t2Residual;
            for (int a = 0; a < NVir; a++)
                for (int b = 0; b < NVir; b++)
                    for (int i = 0; i < NOcc; i++)
                        for (int j = 0; j < NOcc; j++)
                            t2[a, b, i, j] /=
                                mo_occ[i] + mo_occ[j] - mo_vir[a] - mo_vir[b];
        }

        public void ComputeResidual(NDArray eris,
            out NDArray t1Residual, out NDArray t2Residual)
        {
            var eris_oooo = eris[OccRange, OccRange, OccRange, OccRange];
            var eris_oovv = eris[OccRange, OccRange, VirRange, VirRange];
            var eris_vvoo = eris[VirRange, VirRange, OccRange, OccRange];
            var eris_vvvv = eris[VirRange, VirRange, VirRange, VirRange];
            var eris_ovvo = eris[OccRange, VirRange, VirRange, OccRange];
            var eris_ovvv = eris[OccRange, VirRange, VirRange, VirRange];
            var eris_ooov = eris[OccRange, OccRange, OccRange, VirRange];
            var eris_vvvo = eris[VirRange, VirRange, VirRange, OccRange];
            var eris_ovoo = eris[OccRange, VirRange, OccRange, OccRange];

            var tmp = NDArray.ZeroLike(T2);
            Pab(tmp, "ai,bj->abij", T1, T1);
            var ft_vv = NDArray.Contract("fm,mafb->ab", T1, eris_ovvv);
            NDArray.Contract(ft_vv, "afmn,mnbf->ab", -0.5, T2, eris_oovv);
            NDArray.Contract(ft_vv, "afmn,mnbf->ab", -0.25, tmp, eris_oovv);
            var ft_oo = NDArray.Contract("en,inje->ij", T1, eris_ooov);
            NDArray.Contract(ft_oo, "efjn,inef->ij", 0.5, T2, eris_oovv);
            NDArray.Contract(ft_oo, "efjn,inef->ij", 0.25, tmp, eris_oovv);
            var ft_ov = NDArray.Contract("fn,inaf->ia", T1, eris_oovv);
            var w_oooo = eris_oooo.Clone();
            NDArray.Contract(w_oooo, "efij,mnef->mnij", 0.25, T2, eris_oovv);
            NDArray.Contract(w_oooo, "efij,mnef->mnij", 0.25, tmp, eris_oovv);
            Pij(w_oooo, "ej,mnie->mnij", T1, eris_ooov);
            var w_vvvv = eris_vvvv.Clone();
            NDArray.Contract(w_vvvv, "abmn,mnef->abef", 0.25, T2, eris_oovv);
            NDArray.Contract(w_vvvv, "abmn,mnef->abef", 0.25, tmp, eris_oovv);
            Pab(w_vvvv, "bm,maef->abef", T1, eris_ovvv);
            var w_ovvo = eris_ovvo.Clone();
            NDArray.Contract(w_ovvo, "fj,mbef->mbej", T1, eris_ovvv);
            NDArray.Contract(w_ovvo, "bn,mnje->mbej", T1, eris_ooov);
            NDArray.Contract(w_ovvo, "fbjn,mnef->mbej", -0.5, T2, eris_oovv);
            NDArray.Einsum(w_ovvo, "fj,bn,mnef->mbej", -1.0, T1, T1, eris_oovv);

            t1Residual = NDArray.ZeroLike(T1);
            NDArray.Contract(t1Residual, "ei,ae->ai", T1, ft_vv);
            NDArray.Contract(t1Residual, "am,mi->ai", -1.0, T1, ft_oo);
            NDArray.Contract(t1Residual, "aeim,me->ai", T2, ft_ov);
            NDArray.Contract(t1Residual, "fn,nafi->ai", T1, eris_ovvo);
            NDArray.Contract(t1Residual, "efim,maef->ai", -0.5, T2, eris_ovvv);
            NDArray.Contract(t1Residual, "aemn,mnie->ai", -0.5, T2, eris_ooov);
            NDArray.Contract(ft_vv, "bm,me->be", -0.5, T1, ft_ov);
            NDArray.Contract(ft_oo, "ej,me->mj", 0.5, T1, ft_ov);

            t2Residual = eris_vvoo.Clone();
            Pab(t2Residual, "aeij,be->abij", T2, ft_vv);
            Pij(t2Residual, "abim,mj->abij", -1.0, T2, ft_oo);
            NDArray.Contract(t2Residual, "abmn,mnij->abij", 0.5, T2, w_oooo);
            NDArray.Contract(t2Residual, "abmn,mnij->abij", 0.5, tmp, w_oooo);
            NDArray.Contract(t2Residual, "efij,abef->abij", 0.5, T2, w_vvvv);
            NDArray.Contract(t2Residual, "efij,abef->abij", 0.5, tmp, w_vvvv);
            Pabij(t2Residual, "ei,am,mbej->abij", -1.0, T1, T1, eris_ovvo);
            Pabij(t2Residual, "aeim,mbej->abij", T2, w_ovvo);
            Pij(t2Residual, "ei,abej->abij", T1, eris_vvvo);
            Pab(t2Residual, "am,mbij->abij", -1.0, T1, eris_ovoo);
        }

        public void ComputeResidualMin(NDArray eris,
            out NDArray t1Residual, out NDArray t2Residual)
        {
            var eris_oooo = eris[OccRange, OccRange, OccRange, OccRange];
            var eris_ooov = eris[OccRange, OccRange, OccRange, VirRange];
            var eris_oovv = eris[OccRange, OccRange, VirRange, VirRange];
            var eris_ovvv = eris[OccRange, VirRange, VirRange, VirRange];
            var eris_ovvo = eris[OccRange, VirRange, VirRange, OccRange];
            var eris_ovoo = eris[OccRange, VirRange, OccRange, OccRange];
            var eris_vvoo = eris[VirRange, VirRange, OccRange, OccRange];
            var eris_vvvo = eris[VirRange, VirRange, VirRange, OccRange];
            var eris_vvvv = eris[VirRange, VirRange, VirRange, VirRange];

            var nmax = Math.Max(NOcc, NVir);
            tempArray ??= new double[NOcc * NOcc * NVir * NVir];
            var temp_vvoo = new NDArray(tempArray, T2.Lengths);
            NDArray.Contract(0, temp_vvoo, "ai,bj->abij", T1, T1);
            var tmp = temp_vvoo.Clone();
            tmp.SubtractedBy(temp_vvoo.Transpose(1, 0, 2, 3));
            var ft_vv = NDArray.Contract("fm,mafb->ab", T1, eris_ovvv);
            NDArray.Contract(ft_vv, "afmn,mnbf->ab", -0.5, T2, eris_oovv);
            NDArray.Contract(ft_vv, "afmn,mnbf->ab", -0.25, tmp, eris_oovv);
            var ft_oo = NDArray.Contract("en,inje->ij", T1, eris_ooov);
            NDArray.Contract(ft_oo, "efjn,inef->ij", 0.5, T2, eris_oovv);
            NDArray.Contract(ft_oo, "efjn,inef->ij", 0.25, tmp, eris_oovv);
            var ft_ov = NDArray.Contract("fn,inaf->ia", T1, eris_oovv);

            t1Residual = NDArray.ZeroLike(T1);
            NDArray.Contract(t1Residual, "ei,ae->ai", T1, ft_vv);
            NDArray.Contract(t1Residual, "am,mi->ai", -1.0, T1, ft_oo);
            NDArray.Contract(t1Residual, "aeim,me->ai", T2, ft_ov);
            NDArray.Contract(t1Residual, "fn,nafi->ai", T1, eris_ovvo);
            NDArray.Contract(t1Residual, "efim,maef->ai", -0.5, T2, eris_ovvv);
            NDArray.Contract(t1Residual, "aemn,mnie->ai", -0.5, T2, eris_ooov);
            NDArray.Contract(ft_vv, "bm,me->be", -0.5, T1, ft_ov);
            NDArray.Contract(ft_oo, "ej,me->mj", 0.5, T1, ft_ov);

            t2Residual = eris_vvoo.Clone();
            NDArray.Contract(0, temp_vvoo, "aeij,be->abij", T2, ft_vv);
            t2Residual.AddedBy(temp_vvoo).
                SubtractedBy(temp_vvoo.Transpose(1, 0, 2, 3));
            NDArray.Contract(0, temp_vvoo, "abim,mj->abij", T2, ft_oo);
            t2Residual.SubtractedBy(temp_vvoo).
                AddedBy(temp_vvoo.Transpose(0, 1, 3, 2));

            tempArray3 ??= new double[nmax * nmax * nmax * nmax];
            tempArray2 ??= new double[nmax * nmax * nmax * nmax];
            var w_oooo = new NDArray(tempArray3, eris_oooo.Lengths);
            w_oooo.AssignedBy(eris_oooo);
            NDArray.Contract(w_oooo, "efij,mnef->mnij", 0.25, T2, eris_oovv);
            NDArray.Contract(w_oooo, "efij,mnef->mnij", 0.25, tmp, eris_oovv);
            var temp_oooo = new NDArray(tempArray2, w_oooo.Lengths);
            NDArray.Contract(0, temp_oooo, "ej,mnie->mnij", T1, eris_ooov);
            w_oooo.AddedBy(temp_oooo).
                SubtractedBy(temp_oooo.Transpose(0, 1, 3, 2));
            NDArray.Contract(t2Residual, "abmn,mnij->abij", 0.5, T2, w_oooo);
            NDArray.Contract(t2Residual, "abmn,mnij->abij", 0.5, tmp, w_oooo);

            var temp2_vvvv = new NDArray(tempArray2, eris_vvvv.Lengths);
            var w3_vvvv = new NDArray(tempArray3, eris_vvvv.Lengths);
            w3_vvvv.AssignedBy(eris_vvvv);
            NDArray.Contract(w3_vvvv, "abmn,mnef->abef", 0.25, T2, eris_oovv);
            NDArray.Contract(w3_vvvv, "abmn,mnef->abef", 0.25, tmp, eris_oovv);
            NDArray.Contract(0, temp2_vvvv, "bm,maef->abef", T1, eris_ovvv);
            w3_vvvv.AddedBy(temp2_vvvv).
                SubtractedBy(temp2_vvvv.Transpose(1, 0, 2, 3));
            NDArray.Contract(t2Residual, "efij,abef->abij", 0.5, T2, w3_vvvv);
            NDArray.Contract(t2Residual, "efij,abef->abij", 0.5, tmp, w3_vvvv);

            var temp2_vvvo = new NDArray(tempArray2, [NVir, NVir, NVir, NOcc]);
            NDArray.Contract(0, temp2_vvvo, "am,mbej->abej", T1, eris_ovvo);
            NDArray.Contract(0, temp_vvoo, "ei,abej->abij", T1, temp2_vvvo);
            t2Residual.SubtractedBy(temp_vvoo)
                .AddedBy(temp_vvoo.Transpose(1, 0, 2, 3))
                .AddedBy(temp_vvoo.Transpose(0, 1, 3, 2))
                .SubtractedBy(temp_vvoo.Transpose(1, 0, 3, 2));

            var w_ovvo = new NDArray(tempArray3, eris_ovvo.Lengths);
            w_ovvo.AssignedBy(eris_ovvo);
            NDArray.Contract(w_ovvo, "fj,mbef->mbej", T1, eris_ovvv);
            NDArray.Contract(w_ovvo, "bn,mnje->mbej", T1, eris_ooov);
            NDArray.Contract(w_ovvo, "fbjn,mnef->mbej", -0.5, T2, eris_oovv);
            var temp_ovvv = new NDArray(tempArray2, [NOcc, NVir, NVir, NVir]);
            NDArray.Contract(0, temp_ovvv, "bn,mnef->mbef", T1, eris_oovv);
            NDArray.Contract(w_ovvo, "fj,mbef->mbej", -1.0, T1, temp_ovvv);
            NDArray.Contract(0, temp_vvoo, "aeim,mbej->abij", T2, w_ovvo);

            t2Residual.AddedBy(temp_vvoo)
                .SubtractedBy(temp_vvoo.Transpose(1, 0, 2, 3))
                .SubtractedBy(temp_vvoo.Transpose(0, 1, 3, 2))
                .AddedBy(temp_vvoo.Transpose(1, 0, 3, 2));
            NDArray.Contract(0, temp_vvoo, "ei,abej->abij", T1, eris_vvvo);
            t2Residual.AddedBy(temp_vvoo).SubtractedBy(
                temp_vvoo.Transpose(0, 1, 3, 2));
            NDArray.Contract(0, temp_vvoo, "am,mbij->abij", T1, eris_ovoo);
            t2Residual.SubtractedBy(temp_vvoo).AddedBy(
                temp_vvoo.Transpose(1, 0, 2, 3));
        }

        private static void Pab(NDArray self, string expression,
            double scalar, NDArray left, NDArray right)
        {
            var result = NDArray.Contract(expression, scalar, left, right);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([1, 0, 2, 3]));
        }

        private static void Pab(NDArray self, string expression,
            NDArray left, NDArray right)
        {
            var result = NDArray.Contract(expression, left, right);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([1, 0, 2, 3]));
        }

        private static void Pij(NDArray self, string expression,
            double scalar, NDArray left, NDArray right)
        {
            var result = NDArray.Contract(expression, scalar, left, right);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([0, 1, 3, 2]));
        }

        private static void Pij(NDArray self, string expression,
            NDArray left, NDArray right)
        {
            var result = NDArray.Contract(expression, left, right);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([0, 1, 3, 2]));
        }

        private static void Pabij(NDArray self, string expression,
            double scalar, params ReadOnlySpan<NDArray> inputs)
        {
            var result = NDArray.Einsum(expression, scalar, inputs);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([1, 0, 2, 3]))
                .SubtractedBy(result.Transpose([0, 1, 3, 2]))
                .AddedBy(result.Transpose([1, 0, 3, 2]));
        }

        private static void Pabij(NDArray self, string expression,
            NDArray left, NDArray right)
        {
            var result = NDArray.Contract(expression, left, right);
            self.AddedBy(result)
                .SubtractedBy(result.Transpose([1, 0, 2, 3]))
                .SubtractedBy(result.Transpose([0, 1, 3, 2]))
                .AddedBy(result.Transpose([1, 0, 3, 2]));
        }

        public double Energy(NDArray t1, NDArray t2, NDArray eris_vvoo)
        {
            double e_t2 = 0.0;
            double e_t1 = 0.0;
            var t1Mat = t1.AsMatrix();
            for (int a = 0; a < NVir; a++)
            {
                for (int b = 0; b < NVir; b++)
                {
                    for (int i = 0; i < NOcc; i++)
                    {
                        var ai = t1Mat[a, i];
                        for (int j = 0; j < NOcc; j++)
                        {
                            var abij = eris_vvoo[a, b, i, j];
                            e_t2 += t2[a, b, i, j] * abij;
                            e_t1 += ai * t1[b, j] * abij;
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
            var ajkl = NDArray.Contract
                ("abcd,ai->ibcd", ao, coeff);
            var ijcd = NDArray.Contract
                ("ibcd,bj->ijcd", ajkl, coeff);
            var ijkd = NDArray.Contract
                ("ijcd,ck->ijkd", ijcd, coeff);
            var ijkl = NDArray.Contract
                ("ijkd,dl->ijkl", ijkd, coeff);
            return ijkl;
        }

        public static NDArray AO2MO(NDArray ao, NDArray ai, NDArray bj, NDArray ck, NDArray dl)
        {
            if (ao.Rank != 4)
                throw new ArgumentException("AO integral must be a 4D array.");
            GC.Collect();
            var ajkl = NDArray.Contract
                ("abcd,ai->ibcd", ao, ai);
            var ijcd = NDArray.Contract
                ("ibcd,bj->ijcd", ajkl, bj);
            var ijkd = NDArray.Contract
                ("ijcd,ck->ijkd", ijcd, ck);
            var ijkl = NDArray.Contract
                ("ijkd,dl->ijkl", ijkd, dl);
            return ijkl;
        }

        public NDArray MakeErisIncore(RHF mf)
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

        public CCSD(CCD ccd)
        {
            NOcc = ccd.NOcc;
            NVir = ccd.NVir;
            Console.WriteLine($"Nocc={NOcc}, nvir={NVir}");
            NMO = ccd.NMO;
            MOEnergy = ccd.MOEnergies;
            Eris = ccd.Eris;
            MaxCycle = ccd.MaxCycles;
            Tol = ccd.Tol;
            T1 = NDArray.Create([NVir, NOcc]);
            T2 = ccd.T2;
            Eccsd = ccd.Eccd;
            OccRange = new NRange(0, NOcc);
            VirRange = new NRange(NOcc, NMO);
            //DIIS = new CCSDDIIS(NOcc, NVir, 8);
        }

        public CCSD(RHF mf)
        {
            var nao = mf.Envs.NAO;
            var moEnergyArray = NDArray.CreateUninitialized([nao * 2]);
            var unfold = moEnergyArray.Reshape([nao, 2]);
            var aoEne = mf.Es!.AsNDArray();
            unfold[.., 1] = aoEne;
            unfold[.., 0] = aoEne;
            MOEnergy = moEnergyArray.AsVector();
            NOcc = mf.Atoms.Sum(atom => atom.AtomNumber);
            NMO = nao * 2;
            NVir = NMO - NOcc;
            Console.WriteLine($"Nocc={NOcc}, nvir={NVir}");
            Eris = MakeErisIncore(mf);
            Eris = Eris.Transpose(0, 2, 1, 3) - Eris.Transpose(0, 2, 3, 1);
            MaxCycle = 50;
            Tol = 1e-8;
            OccRange = new(0, NOcc);
            VirRange = new(NOcc, NMO);
            //DIIS = new CCSDDIIS(NOcc, NVir, 8);
            var eris_vvoo = Eris[VirRange, VirRange, OccRange, OccRange];
            InitAmps(eris_vvoo, MOEnergy);
        }
    }

    public class CCSDDIIS : DIISBase
    {
        private nint NOcc { get; }
        private nint NVir { get; }

        private NDArray T1s { get; }

        private NDArray T2s { get; }

        public CCSDDIIS(nint nocc, nint nvir, int maxDiis)
            : base(nocc * nocc * nvir * nvir, maxDiis)
        {
            NOcc = nocc;
            NVir = nvir;
            T1s = NDArray.CreateUninitialized([maxDiis, NVir, NOcc]);
            T2s = NDArray.CreateUninitialized([maxDiis, NVir, NVir, NOcc, NOcc]);
        }

        public void Invoke(NDArray t1, NDArray t2, bool update = true)
        {
            var oldT2 = T2s.SliceFirstDim(currentIndex);
            var coeffs = UpdateAndGetCoeffs(error => 
            {
                var errorT2 = error.AsNDArray()
                .Reshape(NVir, NVir, NOcc, NOcc);
                errorT2.AssignedBy(t2).SubtractedBy(oldT2);
            });
            T1s.SliceFirstDim(currentIndex, t1);
            T2s.SliceFirstDim(currentIndex, t2);
            if (coeffs.Length <= 1 || !update)
                return;
            NDArray.BatchAxpy(coeffs, T1s[..currentSize], t1);
            NDArray.BatchAxpy(coeffs, T2s[..currentSize], t2);
        }
    }
}
