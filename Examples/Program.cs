using System.Diagnostics;

namespace Examples;

internal class Program
{
    static void Main()
    {
        MoleculeCase mol = MoleculeCase.CH3OH;

        // RHF
        RHF rhf = new(mol.Atoms, mol.BasisName);
        var totalEnergy = rhf.Run();
        Debug.Assert(IsClose(mol.RhfEnergy, totalEnergy, 1e-8));

        // MP2
        MP2 mP2 = new(rhf);
        var mp2Corr = mP2.Run();
        var mp2Energy = totalEnergy + mp2Corr;
        Debug.Assert(IsClose(mol.Mp2Energy, mp2Energy, 1e-8));

        CCD ccd = new(rhf);
        ccd.Run();
        var ccdEnergy = totalEnergy + ccd.eccd;
        Debug.Assert(IsClose(mol.CCDEnergy, ccdEnergy, 1e-8));

        CCSD ccsd = new(ccd);
        ccsd.Run();
        var ccsdEnergy = totalEnergy + ccsd.eccd;
        Debug.Assert(IsClose(mol.CCSDEnergy, ccsdEnergy, 1e-8));
    }

    static bool IsClose(double a, double b, double tol = 1e-8)
    {
        return Math.Abs(a - b) < tol;
    }
}
