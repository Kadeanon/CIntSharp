using SimpleHelpers;
using SimpleHelpers.LinearAlg;
using SimpleHelpers.MultiAlg;
using SimpleHelpers.MultiAlg.TensorContract;

namespace Examples
{
    internal class MP2
    {
        public RHF RHF { get; }
        public NDArray MOEri { get; }
        public Matrix Coeffs { get; }

        public Vector MOEnergies { get; }
        public nint NOcc { get; }
        public nint NVir { get; }

        public static NDArray compute_mo_eri(Matrix coeffs, NDArray ao_eri)
        {
            var nmo = coeffs.Rows;
            Console.WriteLine("Computing MO ERI...");
            NDArray temp = NDArray.CreateUninitialized
                ([nmo, nmo, nmo, nmo]);
            NDArray mo_eris = NDArray.CreateUninitialized
                ([nmo, nmo, nmo, nmo]);
            NDArray coeffs_view = coeffs.AsNDArray();
            var temp0_view = temp;
            var temp1_view = mo_eris;
            ContractMethods.SimpleContract("ijkl,ip->pjkl", 1.0, ao_eri, coeffs_view, 0.0, temp0_view);
            ContractMethods.SimpleContract("pjkl,jq->pqkl", 1.0, temp0_view, coeffs_view, 0.0, temp1_view);
            ContractMethods.SimpleContract("pqkl,kr->pqrl", 1.0, temp1_view, coeffs_view, 0.0, temp0_view);
            ContractMethods.SimpleContract("pqrl,ls->pqrs", 1.0, temp0_view, coeffs_view, 0.0, temp1_view);
            Console.WriteLine("MO ERI computed successfully.");
            return mo_eris;
        }

        public MP2(RHF rhf)
        {
            RHF = rhf;
            Coeffs = rhf.C!;
            MOEri = compute_mo_eri(Coeffs, rhf.ERI);
            MOEnergies = rhf.Es!;
            nint nmo = Coeffs.Rows;
            NOcc = rhf.Atoms.Sum(atom => atom.AtomNumber) / 2; ;
            NVir = nmo - NOcc;
            Console.WriteLine($"Initialized MP2 with {NOcc} occupied and {NVir} virtual MOs.");
        }


        public double Run()
        {
            Console.WriteLine("Running MP2 energy calculation...");
            Console.WriteLine($"Number of occupied MOs: {NOcc}, Number of virtual MOs: {NVir}");
            double mp2_same_spin = 0.0;
            double mp2_oppo_spin = 0.0;
            var energies = MOEnergies;
            for (int i = 0; i < NOcc; ++i)
            {
                double ei = energies[i];
                for (nint a = NOcc; a < NOcc + NVir; ++a)
                {
                    double ea = energies[a];
                    for (int j = 0; j < NOcc; ++j)
                    {
                        double ej = energies[j];
                        for (nint b = NOcc; b < NOcc + NVir; ++b)
                        {
                            double eb = energies[b];
                            // (ia|jb) * (jb|ia) / (ea + eb - er - es)
                            double eri_arbs = MOEri[i, a, j, b];
                            double eri_asbr = MOEri[i, b, j, a];
                            double denom = ei + ej - ea - eb;
                            if (denom != 0.0)
                            {
                                mp2_same_spin += (eri_arbs - eri_asbr) * eri_arbs / denom;
                                mp2_oppo_spin += eri_arbs * eri_arbs / denom;
                            }
                            else
                            {
                                Console.WriteLine(
                                    "Warning: Division by zero encountered in MP2 energy calculation.");
                            }
                        }
                    }
                }
            }
            Console.WriteLine("MP2 energy contributions calculated successfully.");
            Console.WriteLine($"MP2 same-spin energy contribution:{mp2_same_spin}");
            Console.WriteLine($"MP2 opposite-spin energy contribution:{mp2_oppo_spin}");
            return mp2_same_spin + mp2_oppo_spin;
        }
    }
}
