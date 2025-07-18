using SimpleHelpers.Indices;
using SimpleHelpers.MultiAlg.TensorContract.BSMTC;
using System.Buffers;
using SimpleHelpers.MultiAlg.TensorContract.EinsumTree;
namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {
        public static int DegreeOfParallelism
        {
            get => degreeOfParallelism;
            set => degreeOfParallelism =
                Math.Clamp(value, 1, Environment.ProcessorCount);
        }


        private static int degreeOfParallelism
            = Environment.ProcessorCount / 2;

        public static ArrayPool<nint> ContractNintPool { get; } = ArrayPool<nint>.Create();

        public static void Einsum(NDArray result, string expression, 
            params ReadOnlySpan<NDArray> inputs)
            => result.AddedBy(Einsum(expression, inputs));

        public static void Einsum(NDArray result, string expression, 
            double alpha, params ReadOnlySpan<NDArray> inputs)
            => result.AddedByScaled(alpha, Einsum(expression, inputs));

        public static NDArray Einsum(string expression,
            double alpha, params ReadOnlySpan<NDArray> inputs)
            => Einsum(expression, inputs).ScaledBy(alpha);

        public static NDArray Einsum(string expression,
            params ReadOnlySpan<NDArray> inputs)
        {
            var state = new EinsumState(expression, inputs);
            var ir = state.Parse();
            var result = ir.Invoke(inputs);
            return result;
        }

        internal static int GetLevel(DoubleIndice[] indicesM, DoubleIndice[] indicesN,
            DoubleIndice[] indicesK, double kernelCost)
        {
            double cost = GetCost(
                indicesM, indicesN, indicesK, kernelCost);
            if (cost < 1.25)
                return 1;
            else if (cost < 2)
                return 2;
            else if (cost < 10)
                return 3;
            else
                return 4;
        }

        internal static double GetCost(DoubleIndice[] indicesM, DoubleIndice[] indicesN,
            DoubleIndice[] indicesK, double kernelCost)
        {
            DoubleIndice indiceM = indicesM.Last();
            DoubleIndice indiceN = indicesN.Last();
            DoubleIndice indiceK = indicesK.Last();
            nint SizeM = indicesM.AsSpan().TotalLength();
            nint SizeN = indicesN.AsSpan().TotalLength();
            nint SizeK = indicesK.AsSpan().TotalLength();

            return GetCost(indiceM, indiceN, indiceK,
                SizeM, SizeN, SizeK, kernelCost);
        }

        internal static double GetCost(DoubleIndice indiceM, DoubleIndice indiceN,
            DoubleIndice indiceK, nint sizeM, nint sizeN, nint sizeK, double kernelCost)
        {
            DoubleKernel kernel = new();
            double cost = 0;
            int mc = kernel.mc;
            int mr = kernel.mr;
            int nc = kernel.nc;
            int nr = kernel.nr;
            int kc = kernel.kc;
            nint strideAM = indiceM.AStride;
            nint strideAK = indiceK.AStride;
            nint strideBN = indiceN.AStride;
            nint strideBK = indiceK.BStride;
            nint strideCM = indiceM.BStride;
            nint strideCN = indiceN.BStride;

            double packA = GetCostPerAccess
                (strideAM, strideAK, mc, kc);
            packA /= sizeN;

            double packB = GetCostPerAccess
                (strideBN, strideBK, nc, kc);
            packB *= GetCycles(sizeM, mc);
            packB /= sizeM;

            double unpackC = GetCostPerAccess
                (strideCM, strideCN, mc, nc,
                prefetchRow: true);
            unpackC *= GetCycles(sizeK, kc);
            unpackC /= sizeK;

            double x = kernelCost;
            double kernelFactor = 1;
            if (sizeK < kc)
            {
                x = (kc * kernelFactor + sizeK)
                    / (kernelFactor + 1) / sizeK;
            }

            double parallelFactorM = 1 +
                (DegreeOfParallelism - 1) * 0.9;
            sizeM = Math.Min(sizeM, mc);
            var iBlockNum = (sizeM - 1) / mr + 1;
            var iBlockBatch = (iBlockNum - 1) / 16 + 1;
            var iCoreNum = Math.Min(iBlockBatch, DegreeOfParallelism);
            parallelFactorM /= iCoreNum;

            double parallelFactorN = 1 +
                (DegreeOfParallelism - 1) * 0.9;
            sizeN = Math.Min(sizeN, nc);
            var jBlockNum = (sizeN - 1) / nr + 1;
            var jBlockBatch = (jBlockNum - 1) / 16 + 1;
            var jCoreNum = Math.Min(jBlockBatch, DegreeOfParallelism);
            parallelFactorN /= jCoreNum;

            packA *= parallelFactorM;
            packB *= parallelFactorN;
            unpackC *= parallelFactorM;
            x *= parallelFactorM;

            cost += packA;
            cost += packB;
            cost += unpackC;
            cost += x;
            cost /= (kernelCost + 0.25);
            return cost;

            static double GetCostPerAccess(nint rowStride, nint colStride,
            nint rowBlock, nint colBlock, bool prefetchRow = false)
            {
                double factor = 1.0;
                double simpleFactor = factor * 3;
                int TLBSize = 4096 / 8;
                int cacheLineSize = 64 / 8;

                double rowFactor = 1.0;
                double totalAccess = rowStride * rowBlock;
                double totalPages = totalAccess / TLBSize;
                totalPages = Math.Min(rowBlock, totalPages);
                rowFactor *= totalPages / rowBlock * 300;
                double totalLines = totalAccess / cacheLineSize;
                totalLines = Math.Min(rowBlock, totalLines);
                rowFactor += totalLines / rowBlock * 80;
                if (!prefetchRow)
                {
                    rowFactor *= 0.5;
                }
                else
                {
                    rowFactor *= 0.8;
                }
                factor += rowFactor;

                double colFactor = 1.0;
                totalAccess = colStride * colBlock;
                totalPages = totalAccess / TLBSize;
                totalPages = Math.Min(colBlock, totalPages);
                colFactor *= totalPages / colBlock * 300;
                totalLines = totalAccess / cacheLineSize;
                totalLines = Math.Min(colBlock, totalLines);
                colFactor += totalLines / colBlock * 80;
                colFactor *= 0.8;
                factor += colFactor;

                return Math.Max(simpleFactor, factor);
            }

            static nint GetCycles(nint length, nint block)
            {
                return ((length - 1) / block) + 1;
            }
        }
    }
}
