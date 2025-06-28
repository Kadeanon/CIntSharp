using SimpleHelpers.Indices;
using SimpleHelpers.MultiAlg.TensorContract.SimpleTC;
using System.Buffers;

namespace SimpleHelpers.MultiAlg.TensorContract
{
    public static partial class ContractMethods
    {
        #region Simple TC
        public static NDArray SimpleContract
            (string expression, NDArray left, NDArray right)
            => SimpleContract(expression, 1.0, left, right);

        public static NDArray SimpleContract
            (string expression, double alpha, NDArray left, NDArray right)
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
                char symbol = resultSymbol[i];
                if (indicesA.TryGetValue(symbol, out var input))
                {
                    nint length = input.Length;
                    lengthsC[i] = length;
                }
                else if (indicesB.TryGetValue(symbol, out input))
                {
                    nint length = input.Length;
                    lengthsC[i] = length;
                }
                else
                {
                    throw new ArgumentException($"Symbol '{symbol}' not found in either tensor.");
                }
            }
            NDArray result = NDArray.CreateUninitialized(lengthsC);
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
            var indicesM = indicesAC.Values.ToArray();
            var indicesN = indicesBC.Values.ToArray();
            var indicesK = indicesAB.Values.ToArray();

            if (indicesM.Length == 0)
                indicesM = [new(1, 1, 1)];
            if (indicesN.Length == 0)
                indicesN = [new(1, 1, 1)];
            if (indicesK.Length == 0)
                indicesK = [new(1, 1, 1)];

            SimpleContract(alpha, left, right, 0.0, result, indicesM, indicesN, indicesK);

            return result;
        }

        public static NDArray SimpleContract
            (string expression, double alpha, NDArray left, NDArray right,
            double beta, NDArray result)
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
            var indicesM = indicesAC.Values.ToArray();
            var indicesN = indicesBC.Values.ToArray();
            var indicesK = indicesAB.Values.ToArray();

            if (indicesM.Length == 0)
                indicesM = [new(1, 1, 1)];
            if (indicesN.Length == 0)
                indicesN = [new(1, 1, 1)];
            if (indicesK.Length == 0)
                indicesK = [new(1, 1, 1)];

            SimpleContract(alpha, left, right, beta, result, indicesM, indicesN, indicesK);

            return result;
        }

        private static void SimpleContract(double alpha, NDArray left, NDArray right, 
            double beta, NDArray result,
            DoubleIndice[] indicesM, DoubleIndice[] indicesN, DoubleIndice[] indicesK)
        {
            var indicesMA = indicesM.Select(x => x.A).ToArray();
            var indicesKA = indicesK.Select(x => x.A).ToArray();
            var indicesKB = indicesK.Select(x => x.B).ToArray();
            var indicesNB = indicesN.Select(x => x.A).ToArray();
            var indicesMC = indicesM.Select(x => x.B).ToArray();
            var indicesNC = indicesN.Select(x => x.B).ToArray();

            var matrixA = new SimpleMatrix(left, indicesMA, indicesKA);
            var matrixB = new SimpleMatrix(right, indicesKB, indicesNB);
            var matrixC = new SimpleMatrix(result, indicesMC, indicesNC);

            var m = matrixA.RowLength;
            var n = matrixB.ColLength;
            var k = matrixA.ColLength;

            ParallelOptions options = new() 
            { MaxDegreeOfParallelism = Environment.ProcessorCount / 2 };

            if (beta == 0.0)
            {
                Parallel.For(0, m, options, i =>
                {
                    for (nint j = 0; j < n; j++)
                    {
                        var sum = 0.0;
                        for (nint q = 0; q < k; q++)
                        {
                            sum += matrixA[(nint)i, q] * matrixB[q, j];
                        }
                        matrixC[(nint)i, j] = alpha * sum;
                    }
                });
            }
            else
            {
                Parallel.For(0, m, options, i =>
                {
                    for (nint j = 0; j < n; j++)
                    {
                        var sum = 0.0;
                        for (nint q = 0; q < k; q++)
                        {
                            sum += matrixA[(nint)i, q] * matrixB[q, j];
                        }
                        var last = matrixC[(nint)i, j];
                        matrixC[(nint)i, j] = beta * last + alpha * sum;
                    }
                });
            }
        }

        #endregion Simple TC

    }
}
