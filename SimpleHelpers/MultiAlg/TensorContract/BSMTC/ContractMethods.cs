using SimpleHelpers.Indices;
using SimpleHelpers.MultiAlg.TensorContract.BSMTC;
using System.Buffers;

namespace SimpleHelpers.MultiAlg.TensorContract
{
    public static partial class ContractMethods
    {
        #region BSMTC

        public static void Contract(string expression, double alpha,
            NDArray left, NDArray right, double beta, NDArray result)
        {
            if (beta != 1.0)
                result.ScaledBy(beta);
            Contract(expression, alpha, left, right, result);
        }

        public static void Contract(string expression,
            NDArray left, NDArray right, double beta, NDArray result)
        {
            if (beta != 1.0)
                result.ScaledBy(beta);
            Contract(expression, 1.0, left, right, result);
        }

        public static void Contract(string expression,
            NDArray left, NDArray right, NDArray result)
            => Contract(expression, 1.0, left, right, result);

        public static void Contract(string expression, double alpha,
            NDArray left, NDArray right, NDArray result)
        {
            if (alpha == 0.0)
                return;
            var symbols = expression.Split("->", StringSplitOptions.TrimEntries);
            if (symbols.Length != 2)
                throw new ArgumentException("Expression must contain exactly one '->' symbol.");
            var resultSymbol = symbols[1];
            symbols = symbols[0].Split(",", StringSplitOptions.TrimEntries);
            if (symbols.Length != 2)
                throw new ArgumentException("Expression must contain exactly two symbols before '->'.");
            var leftSymbol = symbols[0];
            var rightSymbol = symbols[1];

            var indicesA = left.Diagonal<NDArray, double>(leftSymbol);
            var indicesB = right.Diagonal<NDArray, double>(rightSymbol);
            var indicesC = result.Diagonal<NDArray, double>(resultSymbol);
            IndiceUtils.Divide(indicesA, indicesB, indicesC,
                out var indicesAB, out var indicesAC, out var indicesBC, out var indicesABC);
            if (indicesA.Count > 0 ||
                indicesB.Count > 0 ||
                indicesC.Count > 0 ||
                indicesABC.Count > 0)
            {
                throw new NotImplementedException("Advanced contracting is not yet implemented.");
            }
            var indicesM = indicesAC.Select(x => x.Value).ToArray();
            var indicesN = indicesBC.Select(x => x.Value).ToArray();
            var indicesK = indicesAB.Select(x => x.Value).ToArray();

            if (indicesM.Length == 0)
                indicesM = [new(1, 1, 1)];
            if (indicesN.Length == 0)
                indicesN = [new(1, 1, 1)];
            if (indicesK.Length == 0)
                indicesK = [new(1, 1, 1)];

            BlockScatterContract(alpha, left, right, result, indicesM, indicesN, indicesK);
        }

        public static NDArray Contract(string expression, NDArray left, NDArray right)
            => Contract(expression, 1.0, left, right);

        public static NDArray Contract(string expression,
            double alpha, NDArray left, NDArray right)
        {
            var symbols = expression.Split("->", StringSplitOptions.TrimEntries);
            if (symbols.Length != 2)
                throw new ArgumentException("Expression must contain exactly one '->' symbol.");
            var resultSymbol = symbols[1];
            symbols = symbols[0].Split(",", StringSplitOptions.TrimEntries);
            if (symbols.Length != 2)
                throw new ArgumentException("Expression must contain exactly two symbols before '->'.");
            var leftSymbol = symbols[0];
            var rightSymbol = symbols[1];

            var indicesA = left.Diagonal<NDArray, double>(leftSymbol);
            var indicesB = right.Diagonal<NDArray, double>(rightSymbol);
            Span<nint> lengthsC = stackalloc nint[resultSymbol.Length];
            for (int i = resultSymbol.Length - 1; i >= 0; i--)
            {
                nint stride = 1;
                char symbol = resultSymbol[i];
                if (indicesA.TryGetValue(symbol, out var input))
                {
                    nint length = input.Length;
                    lengthsC[i] = length;
                    stride *= length;
                }
                else if (indicesB.TryGetValue(symbol, out input))
                {
                    nint length = input.Length;
                    lengthsC[i] = length;
                    stride *= length;
                }
                else
                {
                    throw new ArgumentException($"Symbol '{symbol}' not found in either tensor.");
                }
            }
            NDArray result = NDArray.Create(lengthsC);
            var indicesC = result.Diagonal<NDArray, double>(resultSymbol);
            indicesA = indicesA.Where(kvp => kvp.Value.Length > 1).ToDictionary();
            indicesB = indicesB.Where(kvp => kvp.Value.Length > 1).ToDictionary();
            indicesC = indicesC.Where(kvp => kvp.Value.Length > 1).ToDictionary();
            IndiceUtils.Divide(indicesA, indicesB, indicesC,
                out var indicesAB, out var indicesAC, out var indicesBC, out var indicesABC);
            if (indicesA.Count > 0 ||
                indicesB.Count > 0 ||
                indicesC.Count > 0 ||
                indicesABC.Count > 0
                )
            {
                throw new NotImplementedException("Advanced contracting is not yet implemented.");
            }

            IndiceUtils.Fold(indicesAC);
            var indicesM = indicesAC.Select(x => x.Value).ToArray();
            if (indicesM.Length == 0)
                indicesM = [new(1, 1, 1)];
            Array.Sort(indicesM, (x, y) => y.AStride.CompareTo(x.AStride));
            IndiceUtils.Fold(indicesBC);
            var indicesN = indicesBC.Select(x => x.Value).ToArray();
            if (indicesN.Length == 0)
                indicesN = [new(1, 1, 1)];
            Array.Sort(indicesN, (x, y) => y.AStride.CompareTo(x.AStride));
            IndiceUtils.Fold(indicesAB);
            var indicesK = indicesAB.Select(x => x.Value).ToArray();
            if (indicesK.Length == 0)
                indicesK = [new(1, 1, 1)];
            Array.Sort(indicesK, (x, y) => y.BStride.CompareTo(x.BStride));


            if (indicesM.Last().BStride == 1)
            {
                (left, right) = (right, left);
                (indicesM, indicesN) = (indicesN, indicesM);
                foreach (ref var indice in indicesK.AsSpan())
                {
                    indice = indice.Swap();
                }
            }
            if (alpha != 0.0)
                BlockScatterContract(alpha, left, right, result, indicesM, indicesN, indicesK);

            return result;
        }

        private static void BlockScatterContract(double alpha, NDArray left, NDArray right, NDArray result,
            DoubleIndice[] indicesM, DoubleIndice[] indicesN, DoubleIndice[] indicesK)
        {
            DoubleKernel kernel = default;
            var indicesMA = indicesM.Select(x => x.A).ToArray();
            var indicesKA = indicesK.Select(x => x.A).ToArray();
            var indicesKB = indicesK.Select(x => x.B).ToArray();
            var indicesNB = indicesN.Select(x => x.A).ToArray();
            var indicesMC = indicesM.Select(x => x.B).ToArray();
            var indicesNC = indicesN.Select(x => x.B).ToArray();
            nint mc = kernel.mc, nc = kernel.nc, kc = kernel.kc;
            int mr = kernel.mr, nr = kernel.nr, kr = kernel.kr;

            using var matrixA = new BlockScatterMatrix(left, indicesMA, mr, indicesKA, kr);
            using var matrixB = new BlockScatterMatrix(right, indicesKB, kr, indicesNB, nr);
            using var matrixC = new BlockScatterMatrix(result, indicesMC, mr, indicesNC, nr);
            //matrixB.Transpose();

            var m = matrixA.rowLength;
            var n = matrixB.colLength;
            var k = matrixA.colLength;

            int mMax = (int)Math.Min(mc, m).Align(mr);
            int nMax = (int)Math.Min(nc, n).Align(nr);
            int kMax = (int)Math.Min(kc, k).Align(kr);
            var blockA = ArrayPool<double>.Shared.Rent(kMax * mMax);
            var blockB = ArrayPool<double>.Shared.Rent(kMax * nMax);

            for (nint i = 0; i < m; i += mc)
            {
                int ic = (int)Math.Min(mc, m - i);
                for (nint q = 0; q < k; q += kc)
                {
                    int qc = (int)Math.Min(kc, k - q);
                    int qc_align = qc.Align(kr);
                    var sourceA = matrixA.Slice(i, ic, q, qc);
                    sourceA.Pack(blockA);
                    for (nint j = 0; j < n; j += nc)
                    {
                        int jc = (int)Math.Min(nc, n - j);
                        var sourceB = matrixB.Slice(q, qc, j, jc, trans: true);
                        sourceB.Pack(blockB);
                        var targetC = matrixC.Slice(i, ic, j, jc);
                        KernelParallel kernelParallel =
                            new(targetC, alpha, blockA, blockB, ic, qc_align, jc, kernel);
                        kernelParallel.Invoke();
                    }
                }
            }

            ArrayPool<double>.Shared.Return(blockA);
            ArrayPool<double>.Shared.Return(blockB);
        }

        private static void BSC_MacroKernel(BSMBlock targetC, double alpha,
            double[] blockA, double[] blockB,
            int ic, int qc, int jc, DoubleKernel kernel)
        {
            KernelParallel kernelParallel =
                new(targetC, alpha, blockA, blockB, ic, qc, jc, kernel);
            kernelParallel.Invoke();
        }
        #endregion
    }
}
