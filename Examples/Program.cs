using Examples.Grads;
using Examples.Methods;
using System.Diagnostics;

namespace Examples;

internal class Program
{
    static void Main()
    {
        MoleculeCase mol = MoleculeCase.CH3OH;
        Test(mol);
    }

    static void Test(MoleculeCase mol)
    {
        // RHF
        RHF rhf = new(mol.Atoms, mol.BasisName);
        var totalEnergy = rhf.Run();
        var diff = totalEnergy - mol.RhfEnergy;
        Console.WriteLine($"RHF energy difference: {diff}");
        Debug.Assert(IsClose(mol.RhfEnergy, totalEnergy));

        // MP2
        MP2 mp2 = new(rhf);
        var mp2Energy = mp2.Run() + rhf.TotalEnergy;
        diff = mp2Energy - mol.Mp2Energy;
        Console.WriteLine($"MP2 energy difference: {diff}");
        Debug.Assert(IsClose(mol.Mp2Energy, mp2Energy));

        // CCD
        CCD ccd = new(rhf);
        if (!ccd.Run())
            throw new Exception("CCD failed to converge.");
        var ccdEnergy = ccd.Eccd + rhf.TotalEnergy;
        diff = ccdEnergy - mol.CCDEnergy;
        Console.WriteLine($"CCD energy difference: {diff}");
        Debug.Assert(IsClose(mol.CCDEnergy, ccdEnergy));

        // CCSD
        CCSD ccsd = new(rhf);
        if (!ccsd.Run())
            throw new Exception("CCD failed to converge.");
        var ccsdEnergy = ccsd.Eccsd + rhf.TotalEnergy;
        diff = ccsdEnergy - mol.CCSDEnergy;
        Console.WriteLine($"CCSD energy difference: {diff}");
        Debug.Assert(IsClose(mol.CCSDEnergy, ccsdEnergy));
    }

    static bool IsClose(double a, double b, double tol = 1e-9)
    {
        return Math.Abs(a - b) < tol;
    }
}
