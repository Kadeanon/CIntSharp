using CintSharp.DataStructures;
using CintSharp.DensityFitting;
using CintSharp.DensityFitting.Intor;
using SimpleHelpers;
using SimpleHelpers.LinearAlg;
using SimpleHelpers.MultiAlg;
using SimpleHelpers.MultiAlg.TensorContract;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Numerics.Tensors;
using System.Text;

namespace Examples.Methods.LowScaling;

public class DFRHF
{
    public Atom[] Atoms { get; }

    public DFEnvs Envs { get; }

    public Matrix Ovlp { get; set; }

    public Matrix HCore { get; set; }

    public Matrix Int2c2e { get; set; }

    public NDArray Inv2c2e { get; set; }

    public NDArray Int3c2e { get; set; }

    public Matrix X { get; set; }

    public Matrix XH { get; set; }

    public Matrix P { get; set; }

    public Matrix G { get; set; }

    public Matrix? C { get; set; }

    public Matrix Fock { get; set; }

    public Vector? Es { get; set; }

    public DFRHFDIIS? DIIS { get; set; }

    public const int MaxCycle = 128;
    public const double RequiredDelEnergy = 1e-10;
    private const double RequiredMaxDelF = 1e-9;
    private const double RequiredAverDelF = 1e-10;

    public double ElectronicEnergy
    {
        get
        {
            int nao = Envs.NAO;
            double d = 0;
            for (int i = 0; i < nao; i++)
            {
                for (int j = 0; j < nao; j++)
                {
                    d += P[j, i] * (HCore[i, j] + Fock[i, j]);
                }
            }
            return d / 2;
        }
    }

    public double TotalEnergy
    {
        get
        {
            double result = ElectronicEnergy;
            for (int i = 0; i < Atoms.Length; i++)
            {
                for (int j = i + 1; j < Atoms.Length; j++)
                {
                    var a = Atoms[i];
                    var b = Atoms[j];
                    result += a.AtomNumber * b.AtomNumber
                        / (a.position - b.position).Length;
                }
            }
            return result;
        }
    }

    public DFRHF(IEnumerable<Atom> atoms, string name, string auxName, bool useDIIS = true)
    {
        Atoms = atoms.ToArray();
        Envs = DFEnvs.Create(atoms, atm => name, atm => auxName);
        InitParams();
        DIIS = useDIIS ? new(this) : null;
    }

    public double Run()
    {
        Conver();
        return TotalEnergy;
    }

    [MemberNotNull(
        nameof(Ovlp),
        nameof(HCore),
        nameof(Int2c2e),
        nameof(Inv2c2e),
        nameof(Int3c2e),
        nameof(X),
        nameof(XH),
        nameof(P),
        nameof(G),
        nameof(Fock))]
    public void InitParams()
    {
        int nao = Envs.NAO;

        var ovlpTensor = Envs.GetOvlp();
        Ovlp = ovlpTensor.AsMatrix();

        var hcoreTensor = Envs.GetHCore();
        HCore = hcoreTensor.AsMatrix();
        Int2c2e = Envs.Get2c2e();
        Int3c2e = Envs.Get3c2e();
        var sevd = Int2c2e.SEvd();
        Inv2c2e = sevd.Power(-1).AsNDArray();
        X = Ovlp.SEvd().ReverseSqrt();
        XH = X.Transpose();
        P = Matrix.Create(nao, nao);
        UpdateFockMatrix();
        Fock = HCore + G;
    }

    [MemberNotNull(nameof(G))]
    public void UpdateFockMatrix()
    {
        var PTensor = P.AsNDArray();
        var GTensor = NDArray.Einsum(
            "sl, uvP, PQ, slQ -> uv",
            PTensor, Int3c2e, Inv2c2e, Int3c2e);
        GTensor -= 0.5 * NDArray.Einsum
            ("vl, uvP, PQ, slQ->us",
            PTensor, Int3c2e, Inv2c2e, Int3c2e);
        G = GTensor.AsMatrix();
        Fock = HCore + G;
    }

    [MemberNotNull(nameof(P))]
    private void UpdateDenseMatrix()
    {
        Debug.Assert(C is not null);
        int occ = Atoms.Sum(atom => atom.AtomNumber) / 2;
        P ??= Matrix.Create(C.Rows, C.Cols);
        var pTensor = P.AsSNTensor();
        int rows = (int)pTensor.Lengths[0];
        int cols = (int)pTensor.Lengths[1];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                var ij = 0.0;
                for (int index = 0; index < occ; index++)
                {
                    ij += 2 * C[i, index] * C[j, index];
                }
                pTensor[i, j] = ij;
            }
        }
    }

    public void Conver()
    {
        int cycle = 0;
        bool isConver;
        do
        {
            cycle++;
            Console.WriteLine($"Cycle {cycle}:");
            double oldEnergy = ElectronicEnergy;
            var oldFockMatrix = Fock;
            ConverInner(cycle);
            double delEnergy = oldEnergy - ElectronicEnergy;
            var delFock = (Fock - oldFockMatrix).PointwiseAbs();
            isConver = ConjConver(delEnergy, delFock);
            Console.WriteLine();
        } while (cycle <= MaxCycle && !isConver);
        Console.WriteLine();
        if (!isConver)
        {
            Console.WriteLine($"Not Converged after Cycle {MaxCycle}!");
        }
        else
        {
            Console.WriteLine($"Converged in Cycle {cycle}! Totol Energy= {TotalEnergy}, Electronic Energy= {ElectronicEnergy}.");
        }

    }

    public void ConverInner(int cycle)
    {
        var F_ = XH * Fock * X;
        var evd = F_.SEvd();
        var C_ = evd.EigenVectors;
        Es = evd.EigenValues;
        C = X * C_;
        OutputEnergys();
        UpdateDenseMatrix();
        UpdateFockMatrix();
        TryAddDIIS();
        Console.WriteLine($"Energy: {TotalEnergy}");
    }

    public void TryAddDIIS()
    {
        if (DIIS is null)
            return;
        Fock = DIIS.Invoke(Fock, P);
    }

    public void OutputEnergys()
    {
        Debug.Assert(Es is not null);
        int occ = Atoms.Sum(atom => atom.AtomNumber) / 2;
        int j = 0;
        Console.WriteLine("elec energy: ");
        StringBuilder occStr = new("Occ: [");
        for (; j < occ; j++)
        {
            occStr.Append($"{Es[j]:F4}");
            if (j < occ - 1) occStr.Append(", ");
        }
        occStr.Append(']');
        Console.WriteLine(occStr.ToString());
        StringBuilder virtStr = new("Virt: [");
        for (; j < Es.Length; j++)
        {
            virtStr.Append($"{Es[j]:F4}");
            if (j < Es.Length - 1) virtStr.Append(", ");
        }
        virtStr.Append(']');
        Console.WriteLine(virtStr.ToString());
    }

    public bool ConjConver(double delEnergy, Matrix delFock)
    {
        static string Conj(double result, double need)
        {
            string symbol = result > need ? ">" : "<";
            if (Math.Abs(result) > 1e-4)
            {
                return $"{result:F4}({symbol}{need})";
            }
            return $"{result:E2}({symbol}{need})";
        }
        double dropEnergy = delEnergy;
        double maxDelF = 0;
        double totalSqDelF = 0;
        for (int i = 0; i < delFock.Cols; i++)
        {
            for (int j = 0; j < delFock.Cols; j++)
            {
                double d = delFock[i, j];
                totalSqDelF += d * d;
                if (d > maxDelF)
                {
                    maxDelF = d;
                }
            }
        }
        double averDelF = totalSqDelF / (P.Rows * P.Cols);
        bool isConver = Math.Abs(dropEnergy) < RequiredDelEnergy && maxDelF < RequiredMaxDelF && averDelF < RequiredAverDelF;

        Console.WriteLine($"Energy Drop: {Conj(dropEnergy, RequiredDelEnergy)}, Max DelF: {Conj(maxDelF, RequiredMaxDelF)}, RMS DelF: {Conj(averDelF, RequiredAverDelF)}");
        if (isConver)
        {
            Console.WriteLine("Converged!");
        }
        else
        {
            Console.WriteLine("Not Converged!");
        }
        return isConver;
    }

    public class DFRHFDIIS(DFRHF rhf, int diisSpace = DFRHFDIIS.maxRHFDiis)
        : DIISBase(rhf.Ovlp.Rows * rhf.Ovlp.Rows, diisSpace)
    {
        const int maxRHFDiis = 6;
        readonly Matrix ovlp = rhf.Ovlp;
        protected readonly NDArray inputs = NDArray.Create(
            [diisSpace, rhf.Ovlp.Rows, rhf.Ovlp.Rows]);

        public Matrix Invoke(Matrix fock, Matrix density)
        {
            var B = fock * density * ovlp - ovlp * density * fock;
            var error = B.Flatten();
            var coeffs = UpdateAndGetCoeffs(error);
            inputs[currentIndex, .., ..] = fock.AsNDArray();
            if (coeffs.Length <= 1)
                return fock;
            var fockSlice = inputs[..currentSize];
            NDArray.BatchAxpy(coeffs, fockSlice, fock.AsNDArray());
            return fock;
        }
    }
}
