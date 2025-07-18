using SimpleHelpers.Indices;
using SimpleHelpers.MultiAlg.TensorContract.BSMTC;
using System.Buffers;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {
        #region BSMTC

        public static void Contract(double beta, NDArray result, 
            string expression, double alpha, NDArray left, NDArray right)
        {
            Contract(result.ScaledBy(beta), expression, alpha, left, right);
        }

        public static void Contract(double beta, NDArray result,
            string expression, NDArray left, NDArray right)
        {
            Contract(result.ScaledBy(beta), expression, 1.0, left, right);
        }

        public static void Contract(NDArray result, string expression,
            NDArray left, NDArray right)
            => Contract(result, expression, 1.0, left, right);

        public static void Contract(NDArray result, string expression, 
            double alpha, NDArray left, NDArray right)
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
            Span<TripleIndice> offsets = stackalloc TripleIndice[
                indicesA.Count + indicesB.Count + indicesC.Count + indicesABC.Count];
            int currentIndex = 0;
            foreach (var indice in indicesABC.Values.OrderBy(ind => -ind.CStride))
            {
                offsets[currentIndex++] = indice;
            }
            foreach (var indice in indicesC.Values.OrderBy(ind => -ind.Stride))
            {
                offsets[currentIndex++] = new TripleIndice(
                    indice.Length, 0, 0, indice.Stride);
            }
            foreach (var indice in indicesA.Values.OrderBy(ind => -ind.Stride))
            {
                offsets[currentIndex++] = new TripleIndice(
                    indice.Length, indice.Stride, 0, 0);
            }
            foreach (var indice in indicesB.Values.OrderBy(ind => -ind.Stride))
            {
                offsets[currentIndex++] = new TripleIndice(
                    indice.Length, 0, indice.Stride, 0);
            }

            IndiceUtils.Fold(indicesAC);
            var indicesM = indicesAC.Select(x => x.Value).ToArray();
            Array.Sort(indicesM, (x, y) => y.BStride.CompareTo(x.BStride));
            IndiceUtils.Fold(indicesBC);
            var indicesN = indicesBC.Select(x => x.Value).ToArray();
            Array.Sort(indicesN, (x, y) => y.BStride.CompareTo(x.BStride));
            IndiceUtils.Fold(indicesAB);
            var indicesK = indicesAB.Select(x => x.Value).ToArray();
            Array.Sort(indicesK, (x, y) => y.AStride.CompareTo(x.AStride));

            if (indicesM.Length == 0)
                indicesM = [new(1, 1, 1)];
            if (indicesN.Length == 0)
                indicesN = [new(1, 1, 1)];
            if (indicesK.Length == 0)
                indicesK = [new(1, 1, 1)];

            if (indicesN.Last().BStride == 1)
            {
                (left, right) = (right, left);
                (indicesM, indicesN) = (indicesN, indicesM);
                foreach (ref var indice in indicesK.AsSpan())
                {
                    indice = indice.Swap();
                }
                foreach (ref var indice in offsets)
                {
                    indice =
                        new(indice.Length, indice.BStride, indice.AStride, indice.CStride);
                }
            }

            BlockScatterContract_Silent
                (offsets, alpha, left, 0, right, 0, result, 0, indicesM, indicesN, indicesK);
        }

        public static NDArray Contract(string expression, NDArray left, 
            NDArray right)
            => Contract(expression, 1.0, left, right);

        public static NDArray Contract(string expression,
            double alpha, NDArray left, NDArray right)
        {
            var symbols = expression.Split("->", 
                StringSplitOptions.TrimEntries);
            if (symbols.Length != 2)
                throw new ArgumentException(
                    "Expression must contain exactly one '->' symbol.");
            var resultSymbol = symbols[1];
            symbols = symbols[0].Split(",", 
                StringSplitOptions.TrimEntries);
            if (symbols.Length != 2)
                throw new ArgumentException(
                    "Expression must contain exactly two symbols before '->'.");
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
                    throw new ArgumentException(
                        $"Symbol '{symbol}' not found in either tensor.");
                }
            }
            NDArray result = NDArray.Create(lengthsC);
            var indicesC = result.Diagonal<NDArray, double>(resultSymbol);
            indicesA = indicesA.Where(kvp => kvp.Value.Length > 1).
                ToDictionary();
            indicesB = indicesB.Where(kvp => kvp.Value.Length > 1).
                ToDictionary();
            indicesC = indicesC.Where(kvp => kvp.Value.Length > 1).
                ToDictionary();
            IndiceUtils.Divide(indicesA, indicesB, indicesC,
                out var indicesAB, out var indicesAC, out var indicesBC, 
                out var indicesABC);
            Span<TripleIndice> offsets = stackalloc TripleIndice[
                indicesA.Count + indicesB.Count + indicesC.Count
                + indicesABC.Count];
            int currentIndex = 0;
            foreach (var indice in indicesABC.Values.OrderBy(
                ind => -ind.CStride))
            {
                offsets[currentIndex++] = indice;
            }
            foreach (var indice in indicesC.Values.OrderBy(
                ind => -ind.Stride))
            {
                offsets[currentIndex++] = new TripleIndice(
                    indice.Length, 0, 0, indice.Stride);
            }
            foreach (var indice in indicesA.Values.OrderBy(
                ind => -ind.Stride))
            {
                offsets[currentIndex++] = new TripleIndice(
                    indice.Length, indice.Stride, 0, 0);
            }
            foreach (var indice in indicesB.Values.OrderBy(
                ind => -ind.Stride))
            {
                offsets[currentIndex++] = new TripleIndice(
                    indice.Length, 0, indice.Stride, 0);
            }

            IndiceUtils.Fold(indicesAC);
            var indicesM = indicesAC.Select(x => x.Value).ToArray();
            Array.Sort(indicesM, (x, y) => y.BStride.CompareTo(x.BStride));
            IndiceUtils.Fold(indicesBC);
            var indicesN = indicesBC.Select(x => x.Value).ToArray();
            Array.Sort(indicesN, (x, y) => y.BStride.CompareTo(x.BStride));
            IndiceUtils.Fold(indicesAB);
            var indicesK = indicesAB.Select(x => x.Value).ToArray();
            Array.Sort(indicesK, (x, y) => y.AStride.CompareTo(x.AStride));

            if (indicesM.Length == 0)
                indicesM = [new DoubleIndice(1, 1, 1)];
            if (indicesN.Length == 0)
                indicesN = [new DoubleIndice(1, 1, 1)];
            if (indicesK.Length == 0)
                indicesK = [new DoubleIndice(1, 1, 1)];

            if (indicesN.Last().BStride == 1)
            {
                (left, right) = (right, left);
                (indicesM, indicesN) = (indicesN, indicesM);
                foreach (ref var indice in indicesK.AsSpan())
                {
                    indice = indice.Swap();
                }
                foreach (ref var indice in offsets)
                {
                    indice =
                        new(indice.Length, 
                        indice.BStride, indice.AStride, indice.CStride);
                }
            }

            BlockScatterContract_Silent
                (offsets, alpha, left, 0, right, 0, result, 0, 
                indicesM, indicesN, indicesK);

            return result;
        }

        private static void BlockScatterContract(double alpha, 
            NDArray left, NDArray right, NDArray result,
            DoubleIndice[] indicesM, DoubleIndice[] indicesN, 
            DoubleIndice[] indicesK)
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

            using var matrixA = new BlockScatterMatrix(left, 0, 
                indicesMA, mr, indicesKA, kr);
            using var matrixB = new BlockScatterMatrix(right, 0, 
                indicesKB, kr, indicesNB, nr);
            using var matrixC = new BlockScatterMatrix(result, 0, 
                indicesMC, mr, indicesNC, nr);
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
                        var sourceB = matrixB.Slice(q, qc, j, jc, 
                            trans: true);
                        sourceB.Pack(blockB);
                        var targetC = matrixC.Slice(i, ic, j, jc);
                        KernelParallel kernelParallel =
                            new(targetC, alpha, blockA, blockB, 
                            ic, qc_align, jc, kernel);
                        kernelParallel.Invoke();
                    }
                }
            }

            ArrayPool<double>.Shared.Return(blockA);
            ArrayPool<double>.Shared.Return(blockB);
        }

        private static void BlockScatterContract_Silent
            (ReadOnlySpan<TripleIndice> batchs, double alpha,
            NDArray left, nint leftOffset,
            NDArray right, nint rightOffset,
            NDArray result, nint resultOffset,
            DoubleIndice[] indicesM, DoubleIndice[] indicesN, 
            DoubleIndice[] indicesK)
        {
            if (batchs.Length > 1)
            {
                var currentBatch = batchs[0];
                batchs = batchs[1..];
                for (int i = 0; i < currentBatch.Length; i++)
                {
                    BlockScatterContract_Silent
                        (batchs, alpha, left, leftOffset, 
                        right, rightOffset,
                        result, resultOffset, 
                        indicesM, indicesN, indicesK);
                    leftOffset += currentBatch.AStride;
                    rightOffset += currentBatch.BStride;
                    resultOffset += currentBatch.CStride;
                }
            }
            else
            {
                DoubleKernel kernel = new();
                var indicesMA = indicesM.Select(x => x.A).ToArray();
                var indicesKA = indicesK.Select(x => x.A).ToArray();
                var indicesKB = indicesK.Select(x => x.B).ToArray();
                var indicesNB = indicesN.Select(x => x.A).ToArray();
                var indicesMC = indicesM.Select(x => x.B).ToArray();
                var indicesNC = indicesN.Select(x => x.B).ToArray();
                nint mc = kernel.mc, nc = kernel.nc, kc = kernel.kc;
                int mr = kernel.mr, nr = kernel.nr, kr = kernel.kr;

                var m = indicesMA.AsSpan().TotalLength();
                var n = indicesNB.AsSpan().TotalLength();
                var k = indicesKA.AsSpan().TotalLength();

                int mMax = (int)Math.Min(mc, m).Align(mr);
                int nMax = (int)Math.Min(nc, n).Align(nr);
                int kMax = (int)Math.Min(kc, k).Align(kr);
                nint qAlign = k.Align(kr);
                nint nAlign = n.Align(nr);
                double[] blockA = ArrayPool<double>.Shared.Rent(
                    kMax * mMax);
                double[] blockB = ArrayPool<double>.Shared.Rent(
                    kMax * nMax);

                using var matrixA = new BlockScatterMatrix(left, 
                    leftOffset, indicesMA, mr, indicesKA, kr);
                using var matrixB = new BlockScatterMatrix(right, 
                    rightOffset, indicesKB, kr, indicesNB, nr);
                using var matrixC = new BlockScatterMatrix(result, 
                    resultOffset, indicesMC, mr, indicesNC, nr);

                var currentBatch = batchs.Length == 1 ?
                    batchs[0] : new(1, 0, 0, 0);
                for (int iBatch = 0; iBatch < currentBatch.Length; iBatch++)
                {
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
                                var sourceB = matrixB.Slice(q, qc, j, jc, 
                                    trans: true);
                                sourceB.Pack(blockB);
                                var targetC = matrixC.Slice(i, ic, j, jc);
                                KernelParallel kernelParallel =
                                    new(targetC, alpha, blockA, blockB, 
                                    ic, qc_align, jc, kernel);
                                kernelParallel.Invoke();
                            }
                        }
                    }
                    matrixA.AddOffset(currentBatch.AStride);
                    matrixB.AddOffset(currentBatch.BStride);
                    matrixC.AddOffset(currentBatch.CStride);
                }

                ArrayPool<double>.Shared.Return(blockA);
                ArrayPool<double>.Shared.Return(blockB);
            }
        }

        #endregion
    }
}
