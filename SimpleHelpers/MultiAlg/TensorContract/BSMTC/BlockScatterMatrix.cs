using CommunityToolkit.HighPerformance;
using SimpleHelpers.Indices;
using System.Diagnostics;
using System.Runtime.CompilerServices;

namespace SimpleHelpers.MultiAlg.TensorContract.BSMTC
{
    internal class BlockScatterMatrix(NDArray data,
        SingleIndice[] rowIndices, int rowBlock,
        SingleIndice[] colIndices, int colBlock) : IDisposable
    {
        readonly NDArray data = data;
        #region Rows
        nint[] rowScatters = BuildScatters(rowIndices, rowBlock);
        public nint rowLength = IndiceUtils.TotalLength(rowIndices);
        public int rowBlock = rowBlock;
        #endregion Rows
        #region Cols
        nint[] colScatters = BuildScatters(colIndices, colBlock);
        public nint colLength = IndiceUtils.TotalLength(colIndices);
        public int colBlock = colBlock;

        #endregion Cols

        internal void Pack(Span<double> buffer, nint iStart, int iLength, nint jStart, int jLength, bool trans = false)
        {
            ref double head = ref GetHeadRef();
            int iLengthAlign = iLength.Align(rowBlock);
            int jLengthAlign = jLength.Align(colBlock);
            Debug.Assert(iStart % rowBlock == 0, "iStart must be a multiple of rowBlock.");
            Debug.Assert(iLengthAlign % rowBlock == 0, "iLengthAlign must be a multiple of rowBlock.");
            Debug.Assert(jStart % colBlock == 0, "jStart must be a multiple of colBlock.");
            Debug.Assert(jLengthAlign % colBlock == 0, "jLengthAlign must be a multiple of colBlock.");
            var iScatterOffset = iStart / rowBlock * (rowBlock + 1);
            var jScatterOffset = jStart / colBlock * (colBlock + 1);
            ref var rowScattersHead = ref rowScatters[iScatterOffset];
            ref var colScattersHead = ref colScatters[jScatterOffset];
            var iBlockNum = (iLength - 1) / rowBlock + 1;
            var jBlockNum = (jLength - 1) / colBlock + 1;
            if (trans)
            {
                for (nint j = 0; j < jBlockNum; j++)
                {
                    var jBlockSize = (int)Math.Min(colBlock, jLength - j * colBlock);
                    ref var jToken = ref Unsafe.Add(ref colScattersHead, j * (colBlock + 1));
                    ref var jHead = ref Unsafe.Add(ref jToken, 1);
                    ref var jBufferHead = ref buffer.DangerousGetReferenceAt(
                                j * iLengthAlign * colBlock);
                    for (nint i = 0; i < iBlockNum; i++)
                    {
                        var iBlockSize = (int)Math.Min(rowBlock, iLength - i * rowBlock);
                        ref var iToken = ref Unsafe.Add(ref rowScattersHead, i * (rowBlock + 1));
                        ref var iHead = ref Unsafe.Add(ref iToken, 1);
                        ref var bufferHead = ref Unsafe.Add(ref jBufferHead,
                            i * rowBlock * colBlock);
                        int ii = 0;
                        for (; ii < iBlockSize; ii++)
                        {
                            nint rowOffset = iToken == 0 ?
                                Unsafe.Add(ref iHead, ii) :
                                iHead + ii * iToken;
                            ref double rowPtr = ref Unsafe.Add(ref head, rowOffset);
                            int jj = 0;
                            for (; jj < jBlockSize; jj++)
                            {
                                nint colOffset = jToken == 0 ?
                                    Unsafe.Add(ref jHead, jj) :
                                    jHead + jj * jToken;
                                bufferHead = Unsafe.Add(ref rowPtr, colOffset);
                                bufferHead = ref Unsafe.Add(ref bufferHead, 1);
                            }
                            for (; jj < colBlock; jj++)
                            {
                                bufferHead = 0.0;
                                bufferHead = ref Unsafe.Add(ref bufferHead, 1);
                            }
                        }
                        for (; ii < rowBlock; ii++)
                        {
                            int jj = 0;
                            for (; jj < colBlock; jj++)
                            {
                                bufferHead = 0.0;
                                bufferHead = ref Unsafe.Add(ref bufferHead, 1);
                            }
                        }
                    }
                }
            }
            else
            {
                for (nint i = 0; i < iBlockNum; i++)
                {
                    var iBlockSize = (int)Math.Min(rowBlock, iLength - i * rowBlock);
                    ref var iToken = ref Unsafe.Add(ref rowScattersHead, i * (rowBlock + 1));
                    ref var iHead = ref Unsafe.Add(ref iToken, 1);
                    ref var iBufferHead = ref buffer.DangerousGetReferenceAt(
                                i * jLengthAlign * rowBlock);
                    for (nint j = 0; j < jBlockNum; j++)
                    {
                        var jBlockSize = (int)Math.Min(colBlock, jLength - j * colBlock);
                        ref var jToken = ref Unsafe.Add(ref colScattersHead, j * (colBlock + 1));
                        ref var jHead = ref Unsafe.Add(ref jToken, 1);
                        ref var bufferHead = ref Unsafe.Add(ref iBufferHead,
                            j * colBlock * rowBlock);
                        int jj = 0;
                        for (; jj < jBlockSize; jj++)
                        {
                            nint colOffset = jToken == 0 ?
                                Unsafe.Add(ref jHead, jj) :
                                jHead + jj * jToken;

                            int ii = 0;
                            for (; ii < iBlockSize; ii++)
                            {
                                nint rowOffset = iToken == 0 ?
                                    Unsafe.Add(ref iHead, ii) :
                                    iHead + ii * iToken;
                                ref double rowPtr = ref Unsafe.Add(ref head, rowOffset);

                                bufferHead = Unsafe.Add(ref rowPtr, colOffset);
                                bufferHead = ref Unsafe.Add(ref bufferHead, 1);
                            }
                            for (; ii < rowBlock; ii++)
                            {
                                bufferHead = 0.0;
                                bufferHead = ref Unsafe.Add(ref bufferHead, 1);
                            }
                        }
                        for (; jj < colBlock; jj++)
                        {
                            int ii = 0;
                            for (; ii < rowBlock; ii++)
                            {
                                bufferHead = 0.0;
                                bufferHead = ref Unsafe.Add(ref bufferHead, 1);
                            }
                        }
                    }
                }
            }
        }

        internal nint GetRowOffset(nint index)
        {
            nint head = index / rowBlock * (rowBlock + 1);
            nint left = index % rowBlock;
            ref nint token = ref rowScatters[head];
            ref nint first = ref Unsafe.Add(ref token, 1);
            nint val;
            if (token == 0)
            {
                val = Unsafe.Add(ref first, left);
            }
            else
            {
                val = first + left * token;
            }
            return val;
        }

        internal nint GetColOffset(nint index)
        {
            nint head = index / colBlock * (colBlock + 1);
            nint left = index % colBlock;
            ref nint token = ref colScatters[head];
            ref nint first = ref Unsafe.Add(ref token, 1);
            nint val;
            if (token == 0)
            {
                val = Unsafe.Add(ref first, left);
            }
            else
            {
                val = first + left * token;
            }
            return val;
        }

        public ref double GetHeadRef()
        {
            Span<nint> head = stackalloc nint[data.Rank];
            return ref data[head];
        }

        public void Transpose()
        {
            (rowBlock, colBlock) = (colBlock, rowBlock);
            (rowScatters, colScatters) = (colScatters, rowScatters);
            (rowLength, colLength) = (colLength, rowLength);
        }

        public ref double Local(nint row, nint col)
        {
            nint rowOffset = GetRowOffset(row);
            nint colOffset = GetColOffset(col);
            ref double head = ref GetHeadRef();
            return ref Unsafe.Add(ref head, rowOffset + colOffset);
        }

        public ref double this[nint row, nint col]
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => ref Local(row, col);
        }

        public static nint[] BuildScatters(SingleIndice[] indices, int blockSize)
        {
            nint size = 1;
            for (int i = 0; i < indices.Length; i++)
            {
                size *= indices[i].Length;
            }
            nint blockCount = (size + blockSize - 1) / blockSize;
            var blockScatter = ContractMethods.ContractNintPool.Rent((int)(blockCount) * (blockSize + 1));
            int indiceCount = indices.Length;
            nint[] lengthBuffer = new nint[indiceCount - 1];
            nint[] strideBuffer = new nint[indiceCount - 1];
            ref nint scaRef = ref blockScatter[0];
            var lastIndice = indices[^1];
            nint baseStride = lastIndice.Stride;
            nint baseLength = lastIndice.Length;
            nint times = size / baseLength;
            if (indiceCount >= 2)
            {

                int right = indiceCount - 2;
                ref readonly var last2ndIndice = ref indices[right];
                nint currentScatter = 0;
                nint lastInBlock = 0;
                for (nint i = 0; i < times; i++)
                {
                    nint j = 0;
                    for (; j < lastInBlock; j++)
                    {
                        scaRef = currentScatter;
                        scaRef = ref Unsafe.Add(ref scaRef, 1);
                        currentScatter += baseStride;
                    }
                    for (; j <= baseLength - blockSize; j += blockSize)
                    {
                        scaRef = baseStride;
                        scaRef = ref Unsafe.Add(ref scaRef, 1);
                        scaRef = currentScatter;
                        scaRef = ref Unsafe.Add(ref scaRef, blockSize);
                        currentScatter += baseStride * blockSize;
                    }
                    if (baseLength != j)
                    {
                        //scaRef = 0;//No need
                        scaRef = ref Unsafe.Add(ref scaRef, 1);
                        lastInBlock = (blockSize - (baseLength - j)) % blockSize;
                        for (; j < baseLength; j++)
                        {
                            scaRef = currentScatter;
                            scaRef = ref Unsafe.Add(ref scaRef, 1);
                            currentScatter += baseStride;
                        }
                    }
                    else
                    {
                        lastInBlock = 0;
                    }

                    bool adding = true;
                    currentScatter = 0;
                    ref readonly var indice = ref last2ndIndice;
                    ref nint lenBufRef = ref lengthBuffer.DangerousGetReferenceAt(right);
                    nint lenRef = indice.Length;
                    ref nint strBufRef = ref strideBuffer.DangerousGetReferenceAt(right);
                    nint strRef = indice.Stride;

                    for (var curIndex = right; curIndex >= 0; curIndex--)
                    {
                        if (adding)
                        {
                            lenBufRef++;
                            if (lenBufRef == lenRef)
                            {
                                lenBufRef = 0;
                                strBufRef = 0;
                            }
                            else
                            {
                                adding = false;
                                strBufRef += strRef;
                            }
                            currentScatter += strBufRef;
                            lenBufRef = ref Unsafe.Subtract(ref lenBufRef, 1);
                            indice = ref Unsafe.Subtract(ref Unsafe.AsRef(in indice), 1);
                            lenRef = indice.Length;
                            strRef = indice.Stride;
                        }
                        else
                        {
                            currentScatter += strBufRef;
                        }
                        strBufRef = ref Unsafe.Subtract(ref strBufRef, 1);
                    }
                }
            }
            else
            {
                nint currentScatter = 0;
                nint j = 0;
                for (; j <= baseLength - blockSize; j += blockSize)
                {
                    scaRef = baseStride;
                    scaRef = ref Unsafe.Add(ref scaRef, 1);
                    scaRef = currentScatter;
                    scaRef = ref Unsafe.Add(ref scaRef, blockSize);
                    currentScatter += baseStride * blockSize;
                }
                if (baseLength != j)
                {
                    scaRef = ref Unsafe.Add(ref scaRef, 1);
                    for (; j < baseLength; j++)
                    {
                        scaRef = currentScatter;
                        scaRef = ref Unsafe.Add(ref scaRef, 1);
                        currentScatter += baseStride;
                    }
                }
            }
            return blockScatter;
        }

        public BSMBlock Slice(nint iStart, int iLength, nint jStart, int jLength, bool trans = false)
        {
            ref double head = ref GetHeadRef();
            int iLengthAlign = iLength.Align(rowBlock);
            nint iTotalStart = iStart + iStart / rowBlock;
            int iTotalLength = iLengthAlign + iLengthAlign / rowBlock;
            int jLengthAlign = jLength.Align(colBlock);
            nint jTotalStart = jStart + jStart / colBlock;
            int jTotalLength = jLengthAlign + jLengthAlign / colBlock;
            var colScatterMemory = colScatters.AsMemory((int)jTotalStart, jTotalLength);
            var rowScatterMemory = rowScatters.AsMemory((int)iTotalStart, iTotalLength);
            if (trans)
                return new(data, colScatterMemory, jLength, colBlock,
                    rowScatterMemory, iLength, rowBlock);
            else
                return new(data, rowScatterMemory, iLength, rowBlock,
                    colScatterMemory, jLength, colBlock);
        }

        public void Dispose()
        {
            if (rowScatters != null)
            {
                ContractMethods.ContractNintPool.Return(rowScatters, clearArray: true);
            }
            if (colScatters != null)
            {
                ContractMethods.ContractNintPool.Return(colScatters, clearArray: true);
            }
        }
    }
}
