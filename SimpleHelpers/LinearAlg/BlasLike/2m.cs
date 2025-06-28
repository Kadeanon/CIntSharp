using NetFabric.Numerics.Tensors.Operators;
using System.Runtime.CompilerServices;

namespace SimpleHelpers.LinearAlg
{
    public static partial class BlasLike
    {
        public static void Gemv(double alpha, Matrix a, Vector x, double beta, Vector y)
        {
            nint rows = a.Rows;
            nint cols = a.Cols;
            if ((x.Length != cols) || (y.Length != rows))
                throw new ArgumentException($"Length of x ({x.Length}) must be equal to number of columns a ({cols}).");

            Details.Gemv_Impl(a.Rows, a.Cols, alpha, ref a[0, 0], a.RowStride, a.ColStride, ref x[0], x.Stride, beta, ref y[0], y.Stride);
        }

        public static partial class Details
        {
            public static void Gemv_Impl(nint rows, nint cols, double alpha, ref double aHead, nint aRowStride, nint aColStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                if (rows == 0 || cols == 0)
                    return;
                else if (alpha == 0.0)
                    ApplyWithV_Impl<MultiplyOperator<double>, double>
                        (ref yHead, beta, new(cols, yStride));
                else
                {
                    using var xBuffer = new BufferDVectorSpan(ref xHead, cols, xStride, shouldCopyBack: false);
                    xStride = 1;
                    xHead = ref xBuffer.bufferHead;
                    using var yBuffer = new BufferDVectorSpan(ref yHead, rows, yStride, shouldCopyBack: true);
                    yStride = 1;
                    yHead = ref yBuffer.bufferHead;
                    if (aColStride == 1)
                        Gemv_Kernel_RowMajor_Vector256(rows, cols, alpha, ref aHead, aRowStride, ref xHead, beta, ref yHead);
                    else if (aRowStride == 1)
                        Gemv_Kernel_ColMajor_Vector256(rows, cols, alpha, ref aHead, aColStride, ref xHead, beta, ref yHead);
                    else
                        Gemv_Kernel(rows, cols, alpha, ref aHead, aRowStride, aColStride, ref xHead, beta, ref yHead);
                }
            }

            public static void Gemv_Kernel_RowMajor_Vector256(nint rows, nint cols, double alpha, ref double aHead, nint aRowStride, ref double xHead, double beta, ref double yHead)
            {
                nint pref = DotxF_Kernel_ColMajor_Vector256_PerferredCount;
                nint i = 0;
                for (; i <= rows - pref; i += pref)
                {
                    DotxF_Kernel_ColMajor_Vector256_Perferred_4p6(cols, alpha, ref aHead, aRowStride, ref xHead, 1, beta, ref yHead, yStride: 1);
                    aHead = ref Unsafe.Add(ref aHead, aRowStride * pref);
                    yHead = ref Unsafe.Add(ref yHead, pref);
                }
                for (; i < rows; i++)
                {
                    DotV_Kernel_Vector256(cols, ref aHead, ref xHead, out var rho);
                    yHead = beta * yHead + alpha * rho;
                    aHead = ref Unsafe.Add(ref aHead, aRowStride);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void Gemv_Kernel_ColMajor_Vector256(nint rows, nint cols, double alpha, ref double aHead, nint colStride, ref double xHead, double beta, ref double yHead)
            {
                nint pref = AxpyF_Kernel_ColMajor_Vector256_PerferredCount;
                ApplyWithV_Impl<MultiplyOperator<double>, double>
                    (ref yHead, beta, new(rows, 1));
                nint i = 0;
                for (; i <= cols - pref; i += pref)
                {
                    AxpyF_Kernel_ColMajor_Vector256_Perferred_8p4
                        (rows, alpha, ref aHead, colStride, ref xHead, xStride: 1, 
                        ref yHead, yStride: 1);
                    aHead = ref Unsafe.Add(ref aHead, colStride * pref);
                    xHead = ref Unsafe.Add(ref xHead, pref);
                }
                for (; i < cols; i++)
                {
                    ApplyToWithV_Kernel_Vector<AxpyOperator, double>
                        (rows, ref aHead, alpha * xHead, ref yHead);
                    aHead = ref Unsafe.Add(ref aHead, colStride);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }

            public static void Gemv_Kernel(nint rows, nint cols, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, double beta, ref double yHead)
            {
                for (nint i = 0; i < rows; i++)
                {
                    ref var aRef = ref aHead;
                    ref var xRef = ref xHead;
                    var ySum = 0.0;
                    for (nint j = 0; j < cols; j++)
                    {
                        ySum += aRef * xRef;
                        aRef = ref Unsafe.Add(ref aRef, colStride);
                        xRef = ref Unsafe.Add(ref xRef, 1);
                    }
                    yHead *= beta;
                    yHead += alpha * ySum;
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }
        }

        public static void Ger(double alpha, Vector x, Vector y, Matrix a)
        {
            nint rows = a.Rows;
            nint cols = a.Cols;
            if ((x.Length != rows) || (y.Length != cols))
            {
                throw new ArgumentException($"Length of x ({x.Length}) must be equal to number of rows a ({rows}).");
            }
            Details.Ger_Impl(rows, cols, alpha, ref x[0], x.Stride, ref y[0], y.Stride, ref a[0, 0], a.RowStride, a.ColStride);
        }

        public static partial class Details
        {
            public static void Ger_Impl(nint rows, nint cols, double alpha, ref double xHead, nint xStride, ref double yHead, nint yStride, ref double aHead, nint rowStride, nint colStride)
            {
                if ((rows == 0) || (cols == 0) || (alpha == 0.0))
                    return;
                else
                {
                    using var xBuffer = new BufferDVectorSpan(ref xHead, cols, xStride, shouldCopyBack: false);
                    xStride = 1;
                    xHead = ref xBuffer.bufferHead;
                    using var yBuffer = new BufferDVectorSpan(ref yHead, rows, yStride, shouldCopyBack: false);
                    yStride = 1;
                    yHead = ref yBuffer.bufferHead;
                    if (rowStride == 1)
                        Ger_Kernel_ColMajor_Vector256(rows, cols, alpha, ref xHead, ref yHead, ref aHead, colStride);
                    else if (colStride == 1)
                        Ger_Kernel_RowMajor_Vector256(rows, cols, alpha, ref xHead, ref yHead, ref aHead, rowStride);
                    else
                        Ger_Kernel(rows, cols, alpha, ref xHead, ref yHead, ref aHead, rowStride, colStride);

                }
            }

            private static void Ger_Kernel_ColMajor_Vector256(nint rows, nint cols, double alpha, ref double xHead, ref double yHead, ref double aHead, nint colStride)
            {
                for (nint i = 0; i < cols; i++)
                {
                    ApplyToWithV_Kernel_Vector<AxpyOperator, double>
                        (rows, ref xHead, alpha * yHead, ref aHead);
                    aHead = ref Unsafe.Add(ref aHead, colStride);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            private static void Ger_Kernel_RowMajor_Vector256(nint rows, nint cols, double alpha, ref double xHead, ref double yHead, ref double aHead, nint rowStride)
            {
                for (nint i = 0; i < rows; i++)
                {
                    ApplyToWithV_Kernel_Vector<AxpyOperator, double>
                        (cols, ref yHead, alpha * xHead, ref aHead);
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }

            public static void Ger_Kernel(nint rows, nint cols, double alpha, ref double xHead, ref double yHead, ref double aHead, nint rowStride, nint colStride)
            {
                for (nint i = 0; i < rows; i++)
                {
                    ApplyToWithV_Kernel<AxpyOperator, double>(ref yHead, 
                        alpha * xHead, ref aHead, new(cols, 1, colStride));
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }
        }

        public static void Symv(bool upper, double alpha, Matrix a, Vector x, double beta, Vector y)
        {
            nint length = a.Rows;
            if (length != a.Cols)
                throw new ArgumentException($"Matrix a must be square. Rows: {length}, Cols: {a.Cols}.");
            if ((x.Length != length) || (y.Length != length))
                throw new ArgumentException($"Length of vector must be equal to the length from one dim of matrix.");

            Details.Symv_Impl(upper, length, alpha, ref a[0, 0], a.RowStride, a.ColStride, ref x[0], x.Stride, beta, ref y[0], y.Stride);
        }

        public static partial class Details
        {
            public static void Symv_Impl(bool upper, nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                if (length == 0)
                    return;

                if (upper)
                {
                    (rowStride, colStride) = (colStride, rowStride);
                }
                if (colStride == 1)
                    Symv_Kernel_LoRow_Vector256(length, alpha, ref aHead, rowStride, ref xHead, xStride, beta, ref yHead, yStride);
                else if (rowStride == 1)
                    Symv_Kernel_LoCol_Vector256(length, alpha, ref aHead, colStride, ref xHead, xStride, beta, ref yHead, yStride);
                else 
                    Symv_Kernel_Low(length, alpha, ref aHead, rowStride, colStride, ref xHead, xStride, beta, ref yHead, yStride);
            }

            public static void Symv_Kernel_LoCol_Vector256(nint length, double alpha, ref double aHead, nint colStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride, shouldCopyBack: false);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride, beta, shouldCopyBack: true);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;
                using var yTemp = new BufferDVectorSpan(ref yHead, length, yStride, 0.0);
                ref var yTempHeader = ref yTemp.bufferHead;

                nint pref = DotxAxpyF_Kernel_ColMajor_Vector256_PerferredCount;
                nint i = 0;
                if (i > pref)
                {
                    for (; i <= length - pref; i += pref)
                    {
                        DotxAxpyF_Kernel_ColMajor_Vector256_Perferred_4p4
                            (i, alpha, 
                            ref aHead, colStride,
                            wHead: ref xHead, wStride: 1,
                            xHead: ref xHead, xStride: 1,
                            beta: 1.0,
                            yHead: ref yHead, yStride: 1,
                            zHead: ref yTempHeader, zStride:1);
                        Symv_Kernel_Low
                            (pref, alpha,
                            ref Unsafe.Add(ref aHead, i), rowStride: 1, colStride,
                            ref xHead, xStride,
                            beta: 1.0,
                            ref yTempHeader, yStride: 1);
                        aHead = ref Unsafe.Add(ref aHead, colStride * pref);
                    }
                }
                if(i < length)
                {
                    nint last = length - i;
                    DotxAxpyF_Kernel
                        (i, last, alpha,
                        ref aHead, aRowStride:1, colStride,
                        wHead: ref xHead, wStride: 1,
                        xHead: ref xHead, xStride: 1,
                        1.0,
                        yHead: ref yHead, yStride: 1,
                        zHead: ref yTempHeader, zStride: 1);
                    Symv_Kernel_Low
                        (last, alpha,
                            ref Unsafe.Add(ref aHead, i), rowStride: 1, colStride,
                        ref xHead, xStride,
                        beta,
                        ref yTempHeader, yStride: 1);
                    aHead = ref Unsafe.Add(ref aHead, colStride * pref);
                }
                ApplyToV_Kernel_Vector<AddOperator<double>>
                    (length, ref yTempHeader, ref yHead);
            }

            public static void Symv_Kernel_LoRow_Vector256(nint length, double alpha, ref double aHead, nint rowStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride, shouldCopyBack: false);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride, beta, shouldCopyBack: true);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;
                using var yTemp = new BufferDVectorSpan(ref yHead, length, yStride, 0.0);
                ref var yTempHeader = ref yTemp.bufferHead;

                nint pref = DotxAxpyF_Kernel_RowMajor_Vector256_PerferredCount;
                nint i = 0;
                if (i > pref)
                {
                    for (; i <= length - pref; i += pref)
                    {
                        DotxAxpyF_Kernel_RowMajor_Vector256_Perferred_4p4
                            (i, alpha,
                            ref aHead, rowStride,
                            wHead: ref xHead, wStride: 1,
                            xHead: ref xHead, xStride: 1,
                            beta: 1.0,
                            yHead: ref yHead, yStride: 1,
                            zHead: ref yTempHeader, zStride: 1);
                        Symv_Kernel_Low
                            (pref, alpha,
                            ref Unsafe.Add(ref aHead, i), rowStride, colStride: 1,
                            ref xHead, xStride,
                            beta: 1.0,
                            ref yTempHeader, yStride: 1);
                        aHead = ref Unsafe.Add(ref aHead, rowStride * pref);
                    }
                }
                if (i < length)
                {
                    nint last = length - i;
                    DotxAxpyF_Kernel
                        (i, last, alpha,
                        ref aHead, rowStride, aColStride: 1,
                        wHead: ref xHead, wStride: 1,
                        xHead: ref xHead, xStride: 1,
                        1.0,
                        yHead: ref yHead, yStride: 1,
                        zHead: ref yTempHeader, zStride: 1);
                    Symv_Kernel_Low
                        (last, alpha,
                            ref Unsafe.Add(ref aHead, i), rowStride, colStride: 1,
                        ref xHead, xStride,
                        beta,
                        ref yTempHeader, yStride: 1);
                    aHead = ref Unsafe.Add(ref aHead, rowStride * pref);
                }
                ApplyToV_Kernel_Vector<AddOperator<double>>
                (length, ref yTempHeader, ref yHead);
            }

            public static void Symv_Kernel_Low(nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                ref double xRefI = ref xHead;
                ref double yRefI = ref yHead;
                for (int i = 0; i < length; i++)
                {
                    ref double xRefJ = ref xHead;
                    ref double yRefJ = ref yHead;
                    ref double aRef = ref aHead;

                    yRefI *= beta;
                    for (int j = 0; j < i; j++)
                    {
                        var aVal = alpha * aRef;
                        yRefI += aVal * xRefJ;
                        yRefJ += aVal * xRefI;
                        xRefJ = ref Unsafe.Add(ref xRefJ, xStride);
                        yRefJ = ref Unsafe.Add(ref yRefJ, yStride);
                        aRef = ref Unsafe.Add(ref aRef, colStride);
                    }
                    yRefI += alpha * aRef * xRefI;
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    xRefI = ref Unsafe.Add(ref xRefI, xStride);
                    yRefI = ref Unsafe.Add(ref yRefI, yStride);
                }
            }
        }

        public static void Syr(bool upper, double alpha, Vector x, Matrix a)
        {
            nint length = a.Rows;
            if (length != a.Cols)
                throw new ArgumentException($"Matrix a must be square. Rows: {length}, Cols: {a.Cols}.");
            if (x.Length != length)
                throw new ArgumentException($"Length of vector must be equal to the length from one dim of matrix.");

            Details.Syr_Impl(upper, length, alpha, ref a[0, 0], a.RowStride, a.ColStride, ref x[0], x.Stride);
        }

        public static partial class Details
        {
            public static void Syr_Impl(bool upper, nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride)
            {
                if (length == 0)
                    return;
                if (upper)
                {
                    (rowStride, colStride) = (colStride, rowStride);
                }
                if (colStride == 1)
                    Syr_Kernel_LoRow_Vector256(length, alpha, ref aHead, rowStride, ref xHead, xStride);
                else if (rowStride == 1)
                    Syr_Kernel_LoCol_Vector256(length, alpha, ref aHead, colStride, ref xHead, xStride);
                else
                    Syr_Kernel_Low(length, alpha, ref aHead, rowStride, colStride, ref xHead, xStride);
            }

            public static void Syr_Kernel_LoRow_Vector256(nint length, double alpha, ref double aHead, nint rowStride, ref double xHead, nint xStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                ref var xRefI = ref xBuffer.bufferHead;
                for (nint i = 1; i <= length; i++)
                {
                    ApplyToWithV_Kernel_Vector<AxpyOperator, double>
                    (i, ref xHead, alpha * xRefI, ref aHead);
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    xRefI = ref Unsafe.Add(ref xRefI, 1);
                }
            }

            public static void Syr_Kernel_LoCol_Vector256(nint length, double alpha, ref double aHead, nint colStride, ref double xHead, nint xStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                ref var xRefI = ref xBuffer.bufferHead;
                nint diagStride = colStride + 1;
                for (nint i = length; i > 0; i--)
                {
                    ApplyToWithV_Kernel_Vector<AxpyOperator, double>
                        (i, ref xRefI, alpha * xRefI, ref aHead);
                    aHead = ref Unsafe.Add(ref aHead, diagStride);
                    xRefI = ref Unsafe.Add(ref xRefI, 1);
                }
            }

            public static void Syr_Kernel_Low(nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                ref var xRefI = ref xHead;
                for (nint i = 0; i < length; i++)
                {
                    ref var aRefJ = ref aHead;
                    ref var xRefJ = ref xHead;
                    for (nint j = 0; j < i; j++)
                    {
                        aRefJ += alpha * xRefI * xRefJ;
                        aRefJ = ref Unsafe.Add(ref aRefJ, colStride);
                        xRefJ = ref Unsafe.Add(ref xRefJ, xStride);
                    }
                    aRefJ += alpha * xRefI * xRefI;
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    xRefI = ref Unsafe.Add(ref xRefI, xStride);
                }
            }
        }

        public static void Syr2(bool upper, double alpha, Vector x, Vector y, Matrix a)
        {
            nint length = a.Rows;
            if (length != a.Cols)
                throw new ArgumentException($"Matrix a must be square. Rows: {length}, Cols: {a.Cols}.");
            if ((x.Length != length) || y.Length != length)
                throw new ArgumentException($"Length of vector must be equal to the length from one dim of matrix.");

            Details.Syr2_Impl(upper, length, alpha, ref a[0, 0], a.RowStride, a.ColStride, ref x[0], x.Stride, ref y[0], y.Stride);
        }

        public static partial class Details
        {
            public static void Syr2_Impl(bool upper, nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride, ref double yHead, nint yStride)
            {
                if (length == 0)
                    return;
                if (upper)
                {
                    (rowStride, colStride) = (colStride, rowStride);
                }
                if (colStride == 1)
                    Syr2_Kernel_LoRow_Vector256(length, alpha, ref aHead, rowStride, ref xHead, xStride, ref yHead, yStride);
                else if (rowStride == 1)
                    Syr2_Kernel_LoCol_Vector256(length, alpha, ref aHead, colStride, ref xHead, xStride, ref yHead, yStride);
                else
                    Syr2_Kernel_Low(length, alpha, ref aHead, rowStride, colStride, ref xHead, xStride, ref yHead, yStride);
            }

            public static void Syr2_Kernel_LoRow_Vector256(nint length, double alpha, ref double aHead, nint rowStride, ref double xHead, nint xStride, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;
                ref var xRefI = ref xBuffer.bufferHead;
                ref var yRefI = ref yBuffer.bufferHead;
                for (nint i = 1; i <= length; i++)
                {
                    Axpy2V_Kernel_Vector256(i, alpha * yRefI, alpha * xRefI, ref xHead, ref yHead, ref aHead);
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    xRefI = ref Unsafe.Add(ref xRefI, 1);
                    yRefI = ref Unsafe.Add(ref yRefI, 1);
                }
            }

            public static void Syr2_Kernel_LoCol_Vector256(nint length, double alpha, ref double aHead, nint colStride, ref double xHead, nint xStride, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;
                ref var xRefI = ref xBuffer.bufferHead;
                ref var yRefI = ref yBuffer.bufferHead;
                nint diagStride = colStride + 1;
                for (nint i = length; i > 0; i--)
                {
                    Axpy2V_Kernel_Vector256(i, alpha * yRefI, alpha * xRefI, ref xRefI, ref yRefI, ref aHead);
                    aHead = ref Unsafe.Add(ref aHead, diagStride);
                    xRefI = ref Unsafe.Add(ref xRefI, 1);
                    yRefI = ref Unsafe.Add(ref yRefI, 1);
                }
            }

            public static void Syr2_Kernel_Low(nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;
                ref var xRefI = ref xHead;
                ref var yRefI = ref yHead;
                for (nint i = 0; i < length; i++)
                {
                    ref var aRef = ref aHead;
                    ref var xRefJ = ref xHead;
                    ref var yRefJ = ref yHead;
                    double xValI = xRefI * alpha;
                    double yValI = yRefI * alpha;
                    for (nint j = 0; j <= i; j++)
                    {
                        double del = alpha * xRefI * yRefJ;
                        del += alpha * xRefJ * yRefI; 
                        aRef += del;
                        aRef = ref Unsafe.Add(ref aRef, colStride);
                        xRefJ = ref Unsafe.Add(ref xRefJ, xStride);
                        yRefJ = ref Unsafe.Add(ref yRefJ, xStride);
                    }
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    xRefI = ref Unsafe.Add(ref xRefI, xStride);
                    yRefI = ref Unsafe.Add(ref yRefI, xStride);
                }
            }
        }

        public static void Trmv(bool upper, double alpha, Matrix a, Vector x, double beta, Vector y)
        {
            nint length = a.Rows;
            if (length != a.Cols)
                throw new ArgumentException($"Matrix a must be square. Rows: {length}, Cols: {a.Cols}.");
            if ((x.Length != length) || (y.Length != length))
                throw new ArgumentException($"Length of vector must be equal to the length from one dim of matrix.");

            Details.Trmv_Impl(upper, length, alpha, ref a[0, 0], a.RowStride, a.ColStride, ref x[0], x.Stride, beta, ref y[0], y.Stride);
        }

        public static partial class Details
        {
            public static void Trmv_Impl(bool upper, nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                if (length == 0)
                    return;

                if (upper)
                {
                    if (colStride == 1)
                        Trmv_Kernel_UpRow_Vector256(length, alpha, ref aHead, rowStride, ref xHead, xStride, beta, ref yHead, yStride);
                    else if (rowStride == 1)
                        Trmv_Kernel_UpCol_Vector256(length, alpha, ref aHead, colStride, ref xHead, xStride, beta, ref yHead, yStride);
                    else
                        Trmv_Kernel_Upp(length, alpha, ref aHead, rowStride, colStride, ref xHead, xStride, beta, ref yHead, yStride);
                }
                else
                {
                    if (colStride == 1)
                        Trmv_Kernel_LoRow_Vector256(length, alpha, ref aHead, rowStride, ref xHead, xStride, beta, ref yHead, yStride);
                    else if (rowStride == 1)
                        Trmv_Kernel_LoCol_Vector256(length, alpha, ref aHead, colStride, ref xHead, xStride, beta, ref yHead, yStride);
                    else
                        Trmv_Kernel_Low(length, alpha, ref aHead, rowStride, colStride, ref xHead, xStride, beta, ref yHead, yStride);
                }
            }

            public static void Trmv_Kernel_LoRow_Vector256(nint length, double alpha, ref double aHead, nint rowStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride, shouldCopyBack: false);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride, shouldCopyBack: true);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;

                nint pref = DotxF_Kernel_ColMajor_Vector256_PerferredCount;
                nint i = 0;
                nint iterNums = length / pref;
                nint preIter = length - pref * iterNums;
                Trmv_Kernel_Low
                      (preIter, alpha,
                      ref aHead, rowStride, colStride: 1,
                      ref xHead, xStride: 1,
                      beta,
                      ref yHead, yStride: 1);
                aHead = ref Unsafe.Add(ref aHead, rowStride * preIter);
                yHead = ref Unsafe.Add(ref yHead, preIter);
                i += preIter;
                for (; i < length; i += pref)
                {
                    Trmv_Kernel_Low
                        (pref, alpha,
                        ref Unsafe.Add(ref aHead, i), rowStride, colStride: 1,
                        xHead: ref Unsafe.Add(ref xHead, i), xStride,
                        beta,
                        ref yHead, yStride: 1);
                    DotxF_Kernel_ColMajor_Vector256_Perferred_4p6
                        (i, alpha,
                        ref aHead, rowStride,
                        ref xHead, xStride: 1,
                        beta: 1.0,
                        yHead: ref yHead, yStride: 1);
                    aHead = ref Unsafe.Add(ref aHead, rowStride * pref);
                    yHead = ref Unsafe.Add(ref yHead, pref);
                }
            }

            public static void Trmv_Kernel_LoCol_Vector256(nint length, double alpha, ref double aHead, nint colStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride, shouldCopyBack: false);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride, beta, shouldCopyBack: true);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;

                nint pref = AxpyF_Kernel_ColMajor_Vector256_PerferredCount;
                nint i = 0;
                for (; i <= length - pref; i += pref)
                {
                    Trmv_Kernel_Low
                        (pref, alpha,
                        ref aHead, rowStride: 1, colStride,
                        ref xHead, xStride: 1,
                        beta: 1.0,
                        ref yHead, yStride: 1);
                    aHead = ref Unsafe.Add(ref aHead, pref);
                    yHead = ref Unsafe.Add(ref yHead, pref);
                    AxpyF_Kernel_ColMajor_Vector256_Perferred_8p4(
                        length - pref - i, alpha,
                        ref aHead, colStride,
                        ref xHead, xStride: 1,
                        ref yHead, yStride: 1);
                    xHead = ref Unsafe.Add(ref xHead, pref);
                    aHead = ref Unsafe.Add(ref aHead, colStride * pref);
                }
                nint size = length - i;
                if(size > 0)
                {
                    Trmv_Kernel_Low
                        (size, alpha,
                        ref aHead, rowStride: 1, colStride,
                        ref xHead, xStride: 1,
                        beta: 1.0,
                        ref yHead, yStride: 1);
                }
            }

            public static void Trmv_Kernel_Low(nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride, shouldCopyBack: false);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                for (int i = 0; i < length; i++)
                {
                    ref var aRef = ref aHead;
                    ref var xRef = ref xHead;
                    var ySum = 0.0;
                    for (nint j = 0; j <= i; j++)
                    {
                        ySum += aRef * xRef;
                        aRef = ref Unsafe.Add(ref aRef, colStride);
                        xRef = ref Unsafe.Add(ref xRef, 1);
                    }
                    yHead *= beta;
                    yHead += alpha * ySum;
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    yHead = ref Unsafe.Add(ref yHead, yStride);
                }
            }

            public static void Trmv_Kernel_UpRow_Vector256(nint length, double alpha, ref double aHead, nint rowStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride, shouldCopyBack: false);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride, shouldCopyBack: true);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;

                nint pref = DotxF_Kernel_ColMajor_Vector256_PerferredCount;
                nint i = 0;
                if (length >= pref)
                {
                    for (; i <= length - pref; i += pref)
                    {
                        Trmv_Kernel_Upp
                            (pref, alpha,
                            ref aHead, rowStride, colStride: 1,
                            ref xHead, xStride,
                            beta,
                            ref yHead, yStride: 1);
                        aHead = ref Unsafe.Add(ref aHead, pref);
                        xHead = ref Unsafe.Add(ref xHead, pref);
                        DotxF_Kernel_ColMajor_Vector256_Perferred_4p6
                            (length - i - pref, alpha,
                            ref aHead, rowStride,
                            ref xHead, xStride,
                            beta: 1.0,
                            yHead: ref yHead, yStride: 1);
                        aHead = ref Unsafe.Add(ref aHead, rowStride * pref);
                        yHead = ref Unsafe.Add(ref yHead, pref);
                    }
                }
                if (i < length)
                {
                    nint last = length - i;
                    Trmv_Kernel_Upp
                        (last, alpha,
                        ref aHead, rowStride, colStride: 1,
                        ref xHead, xStride: 1,
                        beta,
                        ref yHead, yStride: 1);
                }
            }

            public static void Trmv_Kernel_UpCol_Vector256(nint length, double alpha, ref double aHead, nint colStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride, shouldCopyBack: false);
                xStride = 1;
                xHead = ref xBuffer.bufferHead;
                using var yBuffer = new BufferDVectorSpan(ref yHead, length, yStride, beta, shouldCopyBack: true);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;

                nint pref = AxpyF_Kernel_ColMajor_Vector256_PerferredCount;
                nint i = 0;
                nint iterNums = length / pref;
                nint preIter = length - pref * iterNums;
                Trmv_Kernel_Upp
                      (preIter, alpha,
                      ref aHead, rowStride: 1, colStride,
                      ref xHead, xStride: 1,
                      beta: 1.0,
                      ref yHead, yStride: 1);
                aHead = ref Unsafe.Add(ref aHead, colStride * preIter);
                xHead = ref Unsafe.Add(ref xHead, preIter);
                i += preIter;
                for (; i < length; i += pref)
                {
                    //AxpyF_Kernel_ColMajor_Vector256_Perferred_8p4
                    //    (i, alpha,
                    //    ref aHead, colStride,
                    //    ref xHead, xStride: 1,
                    //    yHead: ref yHead, yStride: 1);
                    AxpyF_Kernel_ColMajor_Vector256_Perferred_8p4(
                        i, alpha,
                        ref aHead, colStride,
                        ref xHead, xStride: 1,
                        ref yHead, yStride: 1);
                    Trmv_Kernel_Upp
                        (pref, alpha,
                        ref Unsafe.Add(ref aHead, i), rowStride: 1, colStride,
                        ref xHead, xStride,
                        beta: 1.0,
                        ref Unsafe.Add(ref yHead, i), yStride: 1);
                    aHead = ref Unsafe.Add(ref aHead, colStride * pref);
                    xHead = ref Unsafe.Add(ref xHead, pref);
                }
            }

            public static void Trmv_Kernel_Upp(nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double xHead, nint xStride, double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, length, xStride);
                xHead = ref xBuffer.bufferHead;
                nint diagStride = rowStride + colStride;
                for (int i = 0; i < length; i++)
                {
                    ref var aRef = ref aHead;
                    ref var xRef = ref xHead;
                    var ySum = 0.0;
                    for (nint j = i; j < length; j++)
                    {
                        ySum += aRef * xRef;
                        aRef = ref Unsafe.Add(ref aRef, colStride);
                        xRef = ref Unsafe.Add(ref xRef, xStride);
                    }
                    yHead *= beta;
                    yHead += alpha * ySum;
                    aHead = ref Unsafe.Add(ref aHead, diagStride);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, yStride);
                }
            }

        }

        public static void Trsv(bool upper, double alpha, Matrix a, Vector y)
        {
            nint length = a.Rows;
            if (length != a.Cols)
                throw new ArgumentException($"Matrix a must be square. Rows: {length}, Cols: {a.Cols}.");
            if (y.Length != length)
                throw new ArgumentException($"Length of vector must be equal to the length from one dim of matrix.");

            Details.Trsv_Impl(upper, alpha, length, ref a[0, 0], a.RowStride, a.ColStride, ref y[0], y.Stride);
        }

        public static partial class Details
        {
            public static void Trsv_Impl(bool upper, double alpha, nint length, ref double aHead, nint rowStride, nint colStride, ref double yHead, nint yStride)
            {
                if (length == 0)
                    return;

                if (upper)
                {
                    if (colStride == 1)
                        Trsv_Kernel_UpRow_Vector256(length, alpha, ref aHead, rowStride, ref yHead, yStride);
                    else if (rowStride == 1)
                        Trsv_Kernel_UpCol_Vector256(length, alpha, ref aHead, colStride, ref yHead, yStride);
                    else
                        Trsv_Kernel_Upp(length, alpha, ref aHead, rowStride, colStride, ref yHead, yStride);
                }
                else
                {
                    if (colStride == 1)
                        Trsv_Kernel_LoRow_Vector256(length, alpha, ref aHead, rowStride, ref yHead, yStride);
                    else if (rowStride == 1)
                        Trsv_Kernel_LoCol_Vector256(length, alpha, ref aHead, colStride, ref yHead, yStride);
                    else
                        Trsv_Kernel_Low(length, alpha, ref aHead, rowStride, colStride, ref yHead, yStride);
                }
            }

            public static void Trsv_Kernel_LoRow_Vector256(nint length, double alpha, ref double aHead, nint rowStride, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref yHead, length, yStride, shouldCopyBack: false);
                ref var xHead = ref xBuffer.bufferHead;

                nint pref = DotxF_Kernel_ColMajor_Vector256_PerferredCount;
                nint i = 0;
                nint iterNums = length / pref;
                nint preIter = length - pref * iterNums;
                ref var xRef = ref xHead;
                Trsv_Kernel_Low_Inner
                      (preIter, ref aHead, rowStride, colStride: 1, ref xRef);
                aHead = ref Unsafe.Add(ref aHead, rowStride * preIter);
                xRef = ref Unsafe.Add(ref xRef, preIter);
                i += preIter;
                for (; i < length; i += pref)
                {
                    DotxF_Kernel_ColMajor_Vector256_Perferred_4p6
                        (rows:i, alpha: -1.0,
                        ref aHead, rowStride,
                        ref xHead, xStride: 1,
                        beta: 1.0,
                        yHead: ref xRef, yStride: 1);
                    Trsv_Kernel_Low_Inner
                      (pref, ref Unsafe.Add(ref aHead, i), rowStride, colStride: 1, ref xRef);
                    aHead = ref Unsafe.Add(ref aHead, rowStride * pref);
                    xRef = ref Unsafe.Add(ref xRef, pref);
                }
            }

            public static void Trsv_Kernel_LoCol_Vector256(nint length, double alpha, ref double aHead, nint colStride, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref yHead, length, yStride, alpha, shouldCopyBack: true);
                ref var xHead = ref xBuffer.bufferHead;
                ref var xRef = ref xBuffer.bufferHead;
                ref var yRef = ref xBuffer.bufferHead;

                nint pref = AxpyF_Kernel_ColMajor_Vector256_PerferredCount;
                nint i = 0;
                for (; i <= length - pref; i += pref)
                {
                    Trsv_Kernel_Low_Inner
                          (pref, ref aHead, rowStride: 1, colStride, ref xRef);
                    aHead = ref Unsafe.Add(ref aHead, pref);
                    yRef = ref Unsafe.Add(ref yRef, pref);
                    AxpyF_Kernel_ColMajor_Vector256_Perferred_8p4(
                        length - pref - i, alpha: -1.0,
                        ref aHead, colStride,
                        ref xRef, xStride: 1,
                        ref yRef, yStride: 1);
                    aHead = ref Unsafe.Add(ref aHead, colStride * pref);
                    xRef = ref yRef;
                }
                nint size = length - i;
                if (size > 0)
                {
                    Trsv_Kernel_Low_Inner
                          (pref, ref aHead, rowStride: 1, colStride, ref xRef);
                }
            }

            public static void Trsv_Kernel_Low(nint length, double alpha, 
                ref double aHead, nint rowStride, nint colStride, ref double yHead, nint yStride)
            {
                using var xSpan = new BufferDVectorSpan(ref yHead, length, yStride, alpha, shouldCopyBack:true);
                Trsv_Kernel_Low_Inner(length, ref aHead, rowStride, colStride, ref xSpan.bufferHead);
            }

            public static void Trsv_Kernel_Low_Inner(nint length,
                ref double aHead, nint rowStride, nint colStride, ref double xHead)
            {
                ref var xRef = ref xHead;
                ref var aDiagRef = ref aHead;
                nint diagStride = rowStride + colStride;
                for (nint i = 0; i < length; i++)
                {
                    DotV_Impl(ref aHead, ref xHead, new(i, colStride, 1), out var rho);
                    xRef -= rho;
                    xRef /= aDiagRef;
                    aHead = ref Unsafe.Add(ref aHead, rowStride);
                    aDiagRef = ref Unsafe.Add(ref aDiagRef, diagStride);
                    xRef = ref Unsafe.Add(ref xRef, 1);
                }
            }

            public static void Trsv_Kernel_UpRow_Vector256(nint length, double alpha, ref double aHead, nint rowStride, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref yHead, length, yStride, alpha, shouldCopyBack: true);
                ref var xHead = ref xBuffer.bufferHead;

                nint pref = DotxF_Kernel_ColMajor_Vector256_PerferredCount;
                nint iterNums = length / pref;
                nint preIter = length - pref * iterNums;
                ref var xRef = ref Unsafe.Add(ref xHead, pref * iterNums);
                ref var aRef = ref Unsafe.Add(ref aHead, pref * iterNums * (rowStride + 1));
                nint i = 0;
                Trsv_Kernel_Upp_Inner
                      (preIter, ref aRef, rowStride, colStride: 1, ref xRef);
                i += preIter;
                aRef = ref Unsafe.Subtract(ref aRef, rowStride * pref);
                //xRef = ref Unsafe.Subtract(ref xRef, pref);
                ref var yRef = ref Unsafe.Subtract(ref xRef, pref);
                for (; i < length; i += pref)
                {
                    DotxF_Kernel_ColMajor_Vector256_Perferred_4p6
                        (i, alpha: -1.0, ref aRef, rowStride,
                        ref xRef, xStride: 1,
                        beta: 1.0,
                        ref yRef, yStride: 1);
                    aRef = ref Unsafe.Subtract(ref aRef, pref);
                    xRef = ref yRef;
                    Trsv_Kernel_Upp_Inner
                        (pref, ref aRef, rowStride, colStride: 1, ref xRef);
                    aRef = ref Unsafe.Subtract(ref aRef, rowStride * pref);
                    yRef = ref Unsafe.Subtract(ref yRef, pref);
                }
            }

            public static void Trsv_Kernel_UpCol_Vector256(nint length, double alpha, ref double aHead, nint colStride, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref yHead, length, yStride, alpha, shouldCopyBack: true);
                ref var xHead = ref xBuffer.bufferHead;

                nint pref = AxpyF_Kernel_ColMajor_Vector256_PerferredCount;
                nint iterNums = length / pref;
                nint preIter = length - pref * iterNums;
                ref var xRef = ref Unsafe.Add(ref xHead, length - pref);
                ref var aRef = ref Unsafe.Add(ref aHead, (length - pref) * colStride);
                nint i = 0;
                nint mulLength = length;
                for (; i <= length - pref; i += pref)
                {
                    mulLength -= pref;
                    Trsv_Kernel_Upp_Inner(pref, ref Unsafe.Add(ref aRef, mulLength), rowStride: 1, colStride, ref xRef);
                    AxpyF_Kernel_ColMajor_Vector256_Perferred_8p4(
                        mulLength, -1.0,
                        ref aRef, colStride,
                        xHead: ref xRef, xStride: 1,
                        yHead: ref xHead, yStride: 1);
                    xRef = ref Unsafe.Subtract(ref xRef, pref);
                    aRef = ref Unsafe.Subtract(ref aRef, colStride * pref);
                }
                Trsv_Kernel_Upp_Inner(mulLength, ref aHead, rowStride: 1, colStride, ref xHead);
            }

            public static void Trsv_Kernel_Upp(nint length, double alpha, ref double aHead, nint rowStride, nint colStride, ref double yHead, nint yStride)
            {
                using var xSpan = new BufferDVectorSpan(ref yHead, length, yStride, alpha, shouldCopyBack: true);
                Trsv_Kernel_Upp_Inner(length, ref aHead, rowStride, colStride, ref xSpan.bufferHead);
            }

            public static void Trsv_Kernel_Upp_Inner(nint length,
                ref double aHead, nint rowStride, nint colStride, ref double xHead)
            {
                nint diagStride = rowStride + colStride;
                nint index = length - 1;
                ref var aDiagRef = ref Unsafe.Add(ref aHead, index * diagStride);
                ref var xRef = ref Unsafe.Add(ref xHead, index);
                ref var xOldRef = ref Unsafe.Add(ref xRef, 1);
                for (nint i = 0; i < length; i++)
                {
                    DotV_Impl(ref Unsafe.Add(ref aDiagRef, colStride), 
                        yHead: ref xOldRef, new(i, colStride, 1),
                        out var rho);
                    var xRef0 = xRef;
                    xOldRef = ref xRef;
                    xRef -= rho;
                    xRef /= aDiagRef;
                    aDiagRef = ref Unsafe.Subtract(ref aDiagRef, diagStride);
                    xRef = ref Unsafe.Subtract(ref xRef, 1);
                }
            }
        }

    }
}
