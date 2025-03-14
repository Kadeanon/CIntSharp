using CintSharp.DataStructures;
using System.Diagnostics;
using System.Numerics.Tensors;

namespace CintSharp.Examples.GradRHF;

internal class Program
{
    static void Main(string[] args)
    {
        List<Atom> CH3OH = 
        [
            new("C", 0.664617, 0.033231, 0.000000),
            new("H", 1.338596, -1.873145, 0.000000),
            new("H", 1.338631, 0.986406, -1.650963),
            new("H", -1.357391, 0.033256, 0.000000),
            new("O", 1.565402, 1.307100, 2.206427),
            new("H", 0.961143, 3.017646, 2.207215),
        ];
        string BasisName = "sto-3g";
        RHF rhf = new(CH3OH, BasisName);
        var totalEnergy = rhf.Run();
        Debug.Assert(IsClose(-113.54406552303, totalEnergy, 1e-8));
        var grad = new GradRHF(rhf);
        var gradTensor = grad.Run();
        var sum = Tensor.Sum<double>(gradTensor);
        Debug.Assert(IsClose(0.000000, sum, 1e-8));
    }

    static bool IsClose(double a, double b, double tol = 1e-8)
    {
        return Math.Abs(a - b) < tol;
    }
}
