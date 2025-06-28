using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics.X86;

namespace SimpleHelpers.MultiAlg.TensorContract.BSMTC
{
    internal readonly partial struct BSMBlock
    {
        readonly NDArray data;
        #region Rows
        readonly Memory<nint> rowScatters;
        public readonly int rowLength;
        public readonly int rowBlock;
        #endregion Rows
        #region Cols
        readonly Memory<nint> colScatters;
        public readonly int colLength;
        public readonly int colBlock;
        #endregion Cols

        public BSMBlock(NDArray data, Memory<nint> rowScatters, int rowLength, int rowBlock,
            Memory<nint> colScatters, int colLength, int colBlock)
        {
            this.data = data;
            this.rowScatters = rowScatters;
            this.rowLength = rowLength;
            this.rowBlock = rowBlock;
            this.colScatters = colScatters;
            this.colLength = colLength;
            this.colBlock = colBlock;
        }

        public unsafe readonly void Prefetch(int ii, int ir, int jj, int jr, bool perferCol)
        {
            if (!Sse.IsSupported)
            {
                return;
            }
            Span<nint> zeros = stackalloc nint[data.Rank];
            ref double originHead = ref data[zeros];
            var rowSpan = rowScatters.Span;
            var colSpan = colScatters.Span;
            ref nint rowToken = ref rowSpan[ii / rowBlock * (rowBlock + 1)];
            ref nint rowHead = ref Unsafe.Add(ref rowToken, 1);
            ref nint colToken = ref colSpan[jj / colBlock * (colBlock + 1)];
            ref nint colHead = ref Unsafe.Add(ref colToken, 1);
            bool commonRow = rowToken != 0;
            bool commonCol = colToken != 0;

            if (perferCol)
            {
                for (int jjr = 0; jjr < jr; jjr++)
                {
                    var colOffset = colToken == 0 ?
                        Unsafe.Add(ref colHead, jjr) :
                        colHead + colToken * jjr;
                    ref var ptr = ref Unsafe.Add(ref originHead, colOffset);
                    Sse.Prefetch1(Unsafe.AsPointer(
                        ref Unsafe.Add(ref ptr, rowHead)));
                }
            }
            else
            {
                if (commonRow)
                {
                    ref var ptr = ref Unsafe.Add(ref originHead, rowHead + colHead);
                    for (int iir = 0; iir < ir; iir++)
                    {
                        Sse.Prefetch1(Unsafe.AsPointer(
                            ref ptr));
                        ptr = ref Unsafe.Add(ref ptr, rowToken);
                    }
                }
                else
                {
                    ref var ptr = ref Unsafe.Add(ref originHead, colHead);
                    for (int iir = 0; iir < ir; iir++)
                    {
                        Sse.Prefetch1(Unsafe.AsPointer(ref Unsafe.Add(ref originHead,
                            Unsafe.Add(ref rowHead, iir))));
                    }
                }
            }
        }

        public readonly void UnpackScale2(double alpha, int ii, int ir, int jj, int jr, Span<double> buffer, bool perferCol)
        {
            Debug.Assert(buffer.Length >= rowBlock * colBlock);
            Span<nint> zeros = stackalloc nint[data.Rank];
            ref double bufferHead = ref buffer[0];
            ref double originHead = ref data[zeros];
            var rowSpan = rowScatters.Span;
            var colSpan = colScatters.Span;
            ref nint rowToken = ref rowSpan[ii / rowBlock * (rowBlock + 1)];
            ref nint rowHead = ref Unsafe.Add(ref rowToken, 1);
            ref nint colToken = ref colSpan[jj / colBlock * (colBlock + 1)];
            ref nint colHead = ref Unsafe.Add(ref colToken, 1);
            bool commonRow = rowToken != 0;
            bool commonCol = colToken != 0;

            if (perferCol)
            {
                for (int jjr = 0; jjr < jr; jjr++)
                {
                    var colOffset = colToken == 0 ?
                        Unsafe.Add(ref colHead, jjr) :
                        colHead + colToken * jjr;
                    ref var ptr = ref Unsafe.Add(ref originHead, colOffset);
                    ref var bufferPtr = ref bufferHead;
                    for (int iir = 0; iir < ir; iir++)
                    {
                        var rowOffset = rowToken == 0 ?
                            Unsafe.Add(ref rowHead, iir) :
                            rowHead + rowToken * iir;
                        ref var srcPtr = ref Unsafe.Add(ref ptr, rowOffset);
                        srcPtr = alpha * bufferPtr;
                        bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                    }
                    bufferHead = ref Unsafe.Add(ref bufferHead, rowBlock);
                }
            }
            else
            {
                if (commonRow)
                {
                    ref var ptr = ref Unsafe.Add(ref originHead, rowHead);
                    if (commonCol)
                    {
                        for (int iir = 0; iir < ir; iir++)
                        {
                            ref var bufferPtr = ref bufferHead;
                            ref var origPtr = ref Unsafe.Add(ref ptr, colHead);
                            for (int jjr = 0; jjr < jr; jjr++)
                            {
                                origPtr = bufferPtr * alpha;
                                bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                                origPtr = ref Unsafe.Add(ref origPtr, colToken);
                            }
                            bufferHead = ref Unsafe.Add(ref bufferHead, colBlock);
                            ptr = ref Unsafe.Add(ref ptr, rowToken);
                        }
                    }
                    else
                    {
                        for (int iir = 0; iir < ir; iir++)
                        {
                            ref var bufferPtr = ref bufferHead;
                            for (int jjr = 0; jjr < jr; jjr++)
                            {
                                Unsafe.Add(ref ptr, Unsafe.Add(ref colHead, jjr))
                                    = alpha * bufferPtr;
                                bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                            }
                            bufferHead = ref Unsafe.Add(ref bufferHead, colBlock);
                            ptr = ref Unsafe.Add(ref ptr, rowToken);
                        }
                    }
                }
                else
                {
                    if (commonCol)
                    {
                        for (int iir = 0; iir < ir; iir++)
                        {
                            ref var ptr = ref Unsafe.Add(ref originHead,
                                Unsafe.Add(ref rowHead, iir) + colHead);
                            ref var bufferPtr = ref bufferHead;
                            for (int jjr = 0; jjr < jr; jjr++)
                            {
                                ptr = alpha * bufferPtr;
                                bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                                ptr = ref Unsafe.Add(ref ptr, colToken);
                            }
                            bufferHead = ref Unsafe.Add(ref bufferHead, colBlock);
                        }
                    }
                    else
                    {
                        for (int iir = 0; iir < ir; iir++)
                        {
                            ref var ptr = ref Unsafe.Add(ref originHead,
                                Unsafe.Add(ref rowHead, iir));
                            ref var bufferPtr = ref bufferHead;
                            for (int jjr = 0; jjr < jr; jjr++)
                            {
                                Unsafe.Add(ref ptr,
                                    Unsafe.Add(ref colHead, jjr)) =
                                    alpha * bufferPtr;
                                bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                            }
                            bufferHead = ref Unsafe.Add(ref bufferHead, colBlock);
                        }
                    }
                }
            }
        }

        public readonly void UnpackAxpy(double alpha, int ii, int ir, int jj, int jr, Span<double> buffer, bool perferCol)
        {
            Debug.Assert(buffer.Length >= rowBlock * colBlock);
            Span<nint> zeros = stackalloc nint[data.Rank];
            ref double bufferHead = ref buffer[0];
            ref double originHead = ref data[zeros];
            var rowSpan = rowScatters.Span;
            var colSpan = colScatters.Span;
            ref nint rowToken = ref rowSpan[ii / rowBlock * (rowBlock + 1)];
            ref nint rowHead = ref Unsafe.Add(ref rowToken, 1);
            ref nint colToken = ref colSpan[jj / colBlock * (colBlock + 1)];
            ref nint colHead = ref Unsafe.Add(ref colToken, 1);
            bool commonRow = rowToken != 0;
            bool commonCol = colToken != 0;

            if (perferCol)
            {
                for (int jjr = 0; jjr < jr; jjr++)
                {
                    var colOffset = colToken == 0 ?
                        Unsafe.Add(ref colHead, jjr) :
                        colHead + colToken * jjr;
                    ref var ptr = ref Unsafe.Add(ref originHead, colOffset);
                    ref var bufferPtr = ref bufferHead;
                    for (int iir = 0; iir < ir; iir++)
                    {
                        var rowOffset = rowToken == 0 ?
                            Unsafe.Add(ref rowHead, iir) :
                            rowHead + rowToken * iir;
                        ref var srcPtr = ref Unsafe.Add(ref ptr, rowOffset);
                        srcPtr += alpha * bufferPtr;
                        bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                    }
                    bufferHead = ref Unsafe.Add(ref bufferHead, rowBlock);
                }
            }
            else
            {
                if (commonRow)
                {
                    ref var ptr = ref Unsafe.Add(ref originHead, rowHead);
                    if (commonCol)
                    {
                        for (int iir = 0; iir < ir; iir++)
                        {
                            ref var bufferPtr = ref bufferHead;
                            ref var origPtr = ref Unsafe.Add(ref ptr, colHead);
                            for (int jjr = 0; jjr < jr; jjr++)
                            {
                                origPtr += bufferPtr * alpha;
                                bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                                origPtr = ref Unsafe.Add(ref origPtr, colToken);
                            }
                            bufferHead = ref Unsafe.Add(ref bufferHead, colBlock);
                            ptr = ref Unsafe.Add(ref ptr, rowToken);
                        }
                    }
                    else
                    {
                        for (int iir = 0; iir < ir; iir++)
                        {
                            ref var bufferPtr = ref bufferHead;
                            for (int jjr = 0; jjr < jr; jjr++)
                            {
                                Unsafe.Add(ref ptr,Unsafe.Add(ref colHead, jjr))
                                    += alpha * bufferPtr;
                                bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                            }
                            bufferHead = ref Unsafe.Add(ref bufferHead, colBlock);
                            ptr = ref Unsafe.Add(ref ptr, rowToken);
                        }
                    }
                }
                else
                {
                    if (commonCol)
                    {
                        for (int iir = 0; iir < ir; iir++)
                        {
                            ref var ptr = ref Unsafe.Add(ref originHead,
                                Unsafe.Add(ref rowHead, iir) + colHead);
                            ref var bufferPtr = ref bufferHead;
                            for (int jjr = 0; jjr < jr; jjr++)
                            {
                                ptr += alpha * bufferPtr;
                                bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                                ptr = ref Unsafe.Add(ref ptr, colToken);
                            }
                            bufferHead = ref Unsafe.Add(ref bufferHead, colBlock);
                        }
                    }
                    else
                    {
                        for (int iir = 0; iir < ir; iir++)
                        {
                            ref var ptr = ref Unsafe.Add(ref originHead,
                                Unsafe.Add(ref rowHead, iir));
                            ref var bufferPtr = ref bufferHead;
                            for (int jjr = 0; jjr < jr; jjr++)
                            {
                                Unsafe.Add(ref ptr,
                                    Unsafe.Add(ref colHead, jjr)) +=
                                    alpha * bufferPtr;
                                bufferPtr = ref Unsafe.Add(ref bufferPtr, 1);
                            }
                            bufferHead = ref Unsafe.Add(ref bufferHead, colBlock);
                        }
                    }
                }
            }
        }

        internal readonly void Pack(Memory<double> buffer)
        {
            PackParallel pack = new(this, buffer);
            pack.Pack();
        }

        internal readonly nint GetRowOffset(nint index)
        {
            nint head = index / rowBlock * (rowBlock + 1);
            nint left = index % rowBlock;
            ref nint token = ref rowScatters.Span[(int)head];
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
            ref nint token = ref colScatters.Span[(int)head];
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

        public readonly ref double GetHeadRef()
        {
            Span<nint> head = stackalloc nint[data.Rank];
            return ref data[head];
        }

        public readonly ref double Local(nint row, nint col)
        {
            nint rowOffset = GetRowOffset(row);
            nint colOffset = GetColOffset(col);
            ref double head = ref GetHeadRef();
            return ref Unsafe.Add(ref head, rowOffset + colOffset);
        }

        public readonly ref double this[nint row, nint col]
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => ref Local(row, col);
        }
    }
}
