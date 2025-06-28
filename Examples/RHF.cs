using CintSharp.DataStructures;
using CintSharp.Intor;
using SimpleHelpers;
using SimpleHelpers.LinearAlg;
using SimpleHelpers.MultiAlg;
using SimpleHelpers.MultiAlg.TensorContract;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Text;

namespace Examples;

internal class RHF
{
    public Atom[] Atoms { get; }

    public CIntEnvs Envs { get; }

    public Matrix Ovlp { get; set; }

    public Matrix HCore { get; set; }

    public NDArray ERI { get; set; }

    public Matrix X { get; set; }

    public Matrix XH { get; set; }

    public Matrix P { get; set; }

    public Matrix G { get; set; }

    public Matrix? C { get; set; }

    public Matrix Fock { get; set; }

    public Vector? Es { get; set; }

    public RHFDIIS? DIIS { get; set; }

    public const int MaxCycle = 128;
    public const double RequiredDelEnergy = 1e-12;
    private const double RequiredMaxDelF = 1e-12;
    private const double RequiredAverDelF = 1e-12;

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

    public RHF(IEnumerable<Atom> atoms, string name, bool useDIIS = true)
    {
        Atoms = atoms.ToArray();
        Envs = CIntEnvs.Create(atoms, atm => name);
        InitParams();
        DIIS = useDIIS ? new RHFDIIS(this) : null;
    }

    public double Run()
    {
        Conver();
        return TotalEnergy;
    }

    [MemberNotNull(
        nameof(Ovlp),
        nameof(HCore),
        nameof(ERI),
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
        ERI = Envs.GetERI();
        X = Ovlp.SEvd().ReverseSqrt();
        XH = X.Transpose();
        P = Matrix.Create(nao, nao);
        UpdateFockMatrix();
        Fock = HCore + G;
    }

    public static Matrix Tensor2Matrix(NDArray tensor)
    {
        return tensor.AsMatrix();
    }

    [MemberNotNull(nameof(G))]
    public void UpdateFockMatrix()
    {
        nint basisLength = P.Rows;
        G ??= Matrix.Create(basisLength, basisLength);
        var GTensor = G.AsNDArray();
        var PTensor = P.AsNDArray();
        ContractMethods.SimpleContract
            ("sl,uvsl->uv", 1.0, PTensor, ERI, 0.0, GTensor);
        ContractMethods.SimpleContract
            ("ls,mlns->mn", -0.5, PTensor, ERI, 1.0, GTensor);
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
        Fock = DIIS.Update(Fock, P);
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
}
