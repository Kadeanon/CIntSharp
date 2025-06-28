using NetFabric.Numerics.Tensors.Operators;
using SimpleHelpers.Indices;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace SimpleHelpers.LinearAlg
{
    public static partial class BlasLike
    {
        public static void Add(Matrix a, Matrix b)
            => ApplyToM<AddOperator<double>>(a, b);

        public static void Axpy(double alpha, Matrix a, Matrix b)
            => ApplyToWithM<Details.AxpyOperator, double>
            (a, alpha, b);

        public static void Copy(Matrix a, Matrix b)
            => ApplyAssignM<IdentityOperator<double>>(a, b);

        public static void InvScal( double alpha, Matrix a)
        {
            if (alpha == 0.0 || alpha == -0.0 || !double.IsFinite(alpha))
            {
                throw new ArgumentException("Error: alpha must be a finite non-zero value.");
            }
            ApplyWithM<MultiplyOperator<double>, double>(a, 1 / alpha);
        }

        public static void Scal( double alpha, Matrix a)
            => ApplyWithM<MultiplyOperator<double>, double>(a, alpha);

        public static void Scal2( double alpha, Matrix a, Matrix b)
            => ApplyWithAssignM<MultiplyOperator<double>, double>(a, alpha, b);

        public static void Set(double alpha, Matrix a)
            => ApplySelfM<IdentityOperator<double>, double>(a, alpha);

        public static void Sub(Matrix a, Matrix b)
            => ApplyToM<SubtractOperator<double>>(a, b);

        public static void AxpyF( double alpha, Matrix a,  Vector x,  Vector y)
        {
            var rows = a.Rows;
            var cols = a.Cols;
            if (x.Length != cols || y.Length != rows)
                throw new ArgumentException("Matrix and vector dimensions do not match.");
            Details.AxpyF_Impl(rows, cols, alpha, ref a[0, 0], a.RowStride, a.ColStride, ref x[0], x.Stride, ref y[0], y.Stride);
        }

        public static partial class Details
        {
            public static void AxpyF_Impl(nint rows, nint cols,  double alpha, ref double aHead, nint aRowStride, nint aColStride, ref double xHead, nint xStride, ref double yHead, nint yStride)
            {
                if ((rows == 0) || (cols == 0))
                    return;
                else if (cols == 1)
                    ApplyToWithV_Impl<AxpyOperator, double>(ref aHead, 
                        alpha * xHead, ref yHead, new(rows, xStride, yStride));
                else
                {
                    if (aColStride == 1 && rows == AxpyF_Kernel_RowMajor_Vector256_PerferredCount)
                        AxpyF_Kernel_RowMajor_Vector256_Perferred_4p8(rows, alpha, ref aHead, aRowStride, ref xHead, xStride, ref yHead, yStride);
                    else if (aRowStride == 1 && cols == AxpyF_Kernel_ColMajor_Vector256_PerferredCount)
                        AxpyF_Kernel_ColMajor_Vector256_Perferred_8p4(rows, alpha, ref aHead, aColStride, ref xHead, xStride, ref yHead, yStride);
                    else 
                        AxpyF_Kernel(rows, cols, alpha, ref aHead, aRowStride, aColStride, ref xHead, xStride, ref yHead, yStride);
                }
            }

            public const nint AxpyF_Kernel_RowMajor_Vector256_PerferredCount = 8;

            public static void AxpyF_Kernel_RowMajor_Vector256_Perferred_4p8(nint length,  double alpha, ref double aHead, nint aRowStride, ref double xHead, nint xStride, ref double yHead, nint yStride)
            {
                nint iterSize = 4;
                nint i = 0;
                const nint pref = AxpyF_Kernel_RowMajor_Vector256_PerferredCount;
                Span<double> xBuffer = stackalloc double[(int)pref];
                Scal2V_Kernel_Vector256(pref, alpha, ref xHead, ref xBuffer[0]);
                ref var xBuffHead = ref xBuffer[0];
                using var yBuffer = new BufferDVectorSpan(ref yHead, AxpyF_Kernel_RowMajor_Vector256_PerferredCount, yStride, shouldCopyBack: true);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;
                ref var aHead1 = ref Unsafe.Add(ref aHead, 4);

                Vector256<double> xVec0 = Vector256.LoadUnsafe(ref xBuffHead);
                Vector256<double> xVec1 = Vector256.LoadUnsafe(ref xBuffHead, 4);

                if (length >= iterSize)
                {
                    ref var yHead1 = ref Unsafe.Add(ref yHead, 1);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, 2);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, 3);

                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector256<double> aVec00 = Vector256.LoadUnsafe(ref aHead);
                        Vector256<double> aVec01 = Vector256.LoadUnsafe(ref aHead1);
                        aHead = ref Unsafe.Add(ref aHead, aRowStride);
                        aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                        Vector256<double> yVec0 = aVec00 * xVec0;
                        yVec0 = Fma.MultiplyAdd(aVec01, xVec1, yVec0);

                        Vector256<double> aVec10 = Vector256.LoadUnsafe(ref aHead);
                        Vector256<double> aVec11 = Vector256.LoadUnsafe(ref aHead1);
                        aHead = ref Unsafe.Add(ref aHead, aRowStride);
                        aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                        Vector256<double> yVec1 = aVec10 * xVec0;
                        yVec1 = Fma.MultiplyAdd(aVec11, xVec1, yVec0);

                        Vector256<double> aVec20 = Vector256.LoadUnsafe(ref aHead);
                        Vector256<double> aVec21 = Vector256.LoadUnsafe(ref aHead1);
                        aHead = ref Unsafe.Add(ref aHead, aRowStride);
                        aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                        Vector256<double> yVec2 = aVec20 * xVec0;
                        yVec2 = Fma.MultiplyAdd(aVec21, xVec1, yVec0);

                        Vector256<double> aVec30 = Vector256.LoadUnsafe(ref aHead);
                        Vector256<double> aVec31 = Vector256.LoadUnsafe(ref aHead1);
                        aHead = ref Unsafe.Add(ref aHead, aRowStride);
                        aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                        Vector256<double> yVec3 = aVec30 * xVec0;
                        yVec3 = Fma.MultiplyAdd(aVec31, xVec1, yVec0);

                        yHead += yVec0.Sum();
                        yHead1 += yVec1.Sum();
                        yHead2 += yVec2.Sum();
                        yHead3 += yVec3.Sum();

                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                    }
                }

                for (; i < length; i++)
                {
                    Vector256<double> aVec0 = Vector256.LoadUnsafe(ref aHead);
                    Vector256<double> aVec1 = Vector256.LoadUnsafe(ref aHead1);
                    aHead = ref Unsafe.Add(ref aHead, aRowStride);
                    aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                    Vector256<double> yVec = aVec0 * xVec0;
                    yVec = Fma.MultiplyAdd(aVec1, xVec1, yVec);
                    yHead += yVec.Sum();
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }


            public const nint AxpyF_Kernel_ColMajor_Vector256_PerferredCount = 4;

            public static void AxpyF_Kernel_ColMajor_Vector256_Perferred_8p4(nint length,  double alpha, ref double aHead, nint aColStride, ref double xHead, nint xStride, ref double yHead, nint yStride)
            {
                using var yBuffer = new BufferDVectorSpan(ref yHead, AxpyF_Kernel_ColMajor_Vector256_PerferredCount, yStride, shouldCopyBack: true);
                yStride = 1;
                yHead = ref yBuffer.bufferHead;
                nint iterStride = 2;
                nint iterSize = iterStride * Vector256<double>.Count;
                nint i = 0;
                ref double aRef0 = ref aHead;
                ref double aRef1 = ref Unsafe.Add(ref aHead, aColStride);
                ref double aRef2 = ref Unsafe.Add(ref aHead, 2 * aColStride);
                ref double aRef3 = ref Unsafe.Add(ref aHead, 3 * aColStride);
                double fac0 = alpha * xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double fac1 = alpha * xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double fac2 = alpha * xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double fac3 = alpha * xHead;

                if (length >= iterSize)
                {
                    Vector256<double> alpha0 = Vector256.Create(fac0);
                    Vector256<double> alpha1 = Vector256.Create(fac1);
                    Vector256<double> alpha2 = Vector256.Create(fac2);
                    Vector256<double> alpha3 = Vector256.Create(fac3);
                    ref double yHead1 = ref Unsafe.Add(ref yHead, Vector256<double>.Count);

                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector256<double> y0 = Vector256.LoadUnsafe(ref yHead);
                        Vector256<double> y1 = Vector256.LoadUnsafe(ref yHead1);
                        Vector256<double> aVec00 = Vector256.LoadUnsafe(ref aRef0);
                        aRef0 = ref Unsafe.Add(ref aRef0, Vector256<double>.Count);
                        Vector256<double> aVec10 = Vector256.LoadUnsafe(ref aRef0);
                        aRef0 = ref Unsafe.Add(ref aRef0, Vector256<double>.Count);
                        Vector256<double> aVec01 = Vector256.LoadUnsafe(ref aRef1);
                        aRef1 = ref Unsafe.Add(ref aRef1, Vector256<double>.Count);
                        Vector256<double> aVec11 = Vector256.LoadUnsafe(ref aRef1);
                        aRef1 = ref Unsafe.Add(ref aRef1, Vector256<double>.Count);
                        Vector256<double> aVec02 = Vector256.LoadUnsafe(ref aRef2);
                        aRef2 = ref Unsafe.Add(ref aRef2, Vector256<double>.Count);
                        Vector256<double> aVec12 = Vector256.LoadUnsafe(ref aRef2);
                        aRef2 = ref Unsafe.Add(ref aRef2, Vector256<double>.Count);
                        Vector256<double> aVec03 = Vector256.LoadUnsafe(ref aRef3);
                        aRef3 = ref Unsafe.Add(ref aRef3, Vector256<double>.Count);
                        Vector256<double> aVec13 = Vector256.LoadUnsafe(ref aRef3);
                        aRef3 = ref Unsafe.Add(ref aRef3, Vector256<double>.Count);

                        y0 = Fma.MultiplyAdd(aVec00, alpha0, y0);
                        y1 = Fma.MultiplyAdd(aVec10, alpha0, y1);
                        y0 = Fma.MultiplyAdd(aVec01, alpha1, y0);
                        y1 = Fma.MultiplyAdd(aVec11, alpha1, y1);
                        y0 = Fma.MultiplyAdd(aVec02, alpha2, y0);
                        y1 = Fma.MultiplyAdd(aVec12, alpha2, y1);
                        y0 = Fma.MultiplyAdd(aVec03, alpha3, y0);
                        y1 = Fma.MultiplyAdd(aVec13, alpha3, y1);

                        y0.StoreUnsafe(ref yHead);
                        y1.StoreUnsafe(ref yHead1);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                    }
                }

                for (; i < length; i++)
                {
                    double yVal = yHead;
                    yVal += fac0 * aRef0;
                    yVal += fac1 * aRef1;
                    yVal += fac2 * aRef2;
                    yVal += fac3 * aRef3;
                    yHead = yVal;

                    yHead = ref Unsafe.Add(ref yHead, 1);
                    aRef0 = ref Unsafe.Add(ref aRef0, 1);
                    aRef1 = ref Unsafe.Add(ref aRef1, 1);
                    aRef2 = ref Unsafe.Add(ref aRef2, 1);
                    aRef3 = ref Unsafe.Add(ref aRef3, 1);
                }
            }
            
            public static void AxpyF_Kernel(nint rows, nint cols,  double alpha, ref double aHead, nint aRowStride, nint aColStride, ref double xHead, nint xStride, ref double yHead, nint yStride)
            {
                for (nint i = 0; i < cols; i++)
                {
                    ApplyToWithV_Impl<AxpyOperator, double>(ref aHead, alpha * xHead, 
                        ref yHead, new(rows, aRowStride, yStride));
                    aHead = ref Unsafe.Add(ref aHead, aColStride);
                    xHead = ref Unsafe.Add(ref xHead, xStride);
                }
            }

            private static void Scal2V_Kernel_Vector256(nint length, double alpha, ref double xHead, ref double yHead)
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector256<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector256<double>.Count);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, Vector256<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector256<double>.Count * 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, Vector256<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector256<double>.Count * 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, Vector256<double>.Count * 3);
                    Vector256<double> alphaVec = Vector256.Create(alpha);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector256<double> xVec0 = Vector256.LoadUnsafe(ref xHead);
                        Vector256<double> xVec1 = Vector256.LoadUnsafe(ref xHead1);
                        Vector256<double> xVec2 = Vector256.LoadUnsafe(ref xHead2);
                        Vector256<double> xVec3 = Vector256.LoadUnsafe(ref xHead3);
                        xVec0 *= alphaVec;
                        xVec1 *= alphaVec;
                        xVec2 *= alphaVec;
                        xVec3 *= alphaVec;
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                        xVec0.StoreUnsafe(ref yHead);
                        xVec1.StoreUnsafe(ref yHead1);
                        xVec2.StoreUnsafe(ref yHead2);
                        xVec3.StoreUnsafe(ref yHead3);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }
        }

        public static void DotxF( double alpha, Matrix a,  Vector x,  double beta,  Vector y)
        {
            if (alpha == 0)
            {
                Scal(beta, y);
                return;
            }
            var rows = a.Rows;
            var cols = a.Cols;
            if (y.Length != cols || x.Length != rows)
                throw new ArgumentException("Matrix and vector dimensions do not match.");
            Details.DotxF_Impl(rows, cols, alpha, ref a[0, 0], a.RowStride, a.ColStride, ref x[0], x.Stride, beta, ref y[0], y.Stride);
        }
        
        public static partial class Details
        {
            public static void DotxF_Impl(nint rows, nint cols,  double alpha, ref double aHead, nint aRowStride, nint aColStride, ref double xHead, nint xStride,  double beta, ref double yHead, nint yStride)
            {
                if ((rows == 0) || (cols == 0))
                    return;
                else if (cols == 1)
                {
                    yHead *= beta;
                    DotV_Impl(ref aHead, ref xHead, new(rows, aRowStride, xStride), out var rho);
                    yHead += alpha * rho;
                }
                else if (aColStride == 1 && cols == DotxF_Kernel_RowMajor_Vector256_PerferredCount)
                    DotxF_Kernel_RowMajor_Vector256_Perferred_4p8(rows, alpha, ref aHead, aRowStride, ref xHead, xStride, beta, ref yHead, yStride);
                else if(aRowStride == 1 && cols == DotxF_Kernel_ColMajor_Vector256_PerferredCount)
                    DotxF_Kernel_ColMajor_Vector256_Perferred_4p6(rows, alpha, ref aHead, aColStride, ref xHead, xStride, beta, ref yHead, yStride);
                else
                    DotxF_Kernel(rows, cols, alpha, ref aHead, aRowStride, aColStride, ref xHead, xStride, beta, ref yHead, yStride);
            }

            public const nint DotxF_Kernel_RowMajor_Vector256_PerferredCount = 8;

            public static void DotxF_Kernel_RowMajor_Vector256_Perferred_4p8(nint rows,  double alpha, ref double aHead, nint aRowStride, ref double xHead, nint xStride,  double beta, ref double yHead, nint yStride)
            {
                nint iterSize = 4;
                Span<double> yBuffer = stackalloc double[(int)DotxF_Kernel_RowMajor_Vector256_PerferredCount];
                ref var yBufferHead = ref MemoryMarshal.GetReference(yBuffer);
                Vector256<double> yVec0 = Vector256.LoadUnsafe(ref yBufferHead);
                Vector256<double> yVec1 = Vector256.LoadUnsafe(ref yBufferHead, 4);

                ref var aHead1 = ref Unsafe.Add(ref aHead, 4);
                nint i = 0;
                for (; i <= rows - iterSize; i += iterSize)
                {
                    double xScalar0 = xHead * alpha;
                    Vector256<double> xVec0 = Vector256.Create(xScalar0);
                    Vector256<double> aVec00 = Vector256.LoadUnsafe(ref aHead);
                    Vector256<double> aVec10 = Vector256.LoadUnsafe(ref aHead1);
                    xHead = ref Unsafe.Add(ref xHead, xStride);
                    aHead = ref Unsafe.Add(ref aHead, aRowStride);
                    aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                    yVec0 = Fma.MultiplyAdd(aVec00, xVec0, yVec0);
                    yVec1 = Fma.MultiplyAdd(aVec10, xVec0, yVec1);

                    double xScalar1 = xHead * alpha;
                    Vector256<double> xVec1 = Vector256.Create(xScalar1);
                    Vector256<double> aVec01 = Vector256.LoadUnsafe(ref aHead);
                    Vector256<double> aVec11 = Vector256.LoadUnsafe(ref aHead1);
                    xHead = ref Unsafe.Add(ref xHead, xStride);
                    aHead = ref Unsafe.Add(ref aHead, aRowStride);
                    aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                    yVec0 = Fma.MultiplyAdd(aVec01, xVec1, yVec0);
                    yVec1 = Fma.MultiplyAdd(aVec11, xVec1, yVec1);

                    double xScalar2 = xHead * alpha;
                    Vector256<double> xVec2 = Vector256.Create(xScalar2);
                    Vector256<double> aVec02 = Vector256.LoadUnsafe(ref aHead);
                    Vector256<double> aVec12 = Vector256.LoadUnsafe(ref aHead1);
                    xHead = ref Unsafe.Add(ref xHead, xStride);
                    aHead = ref Unsafe.Add(ref aHead, aRowStride);
                    aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                    yVec0 = Fma.MultiplyAdd(aVec02, xVec2, yVec0);
                    yVec1 = Fma.MultiplyAdd(aVec12, xVec2, yVec1);
                    
                    double xScalar3 = xHead * alpha;
                    Vector256<double> xVec3 = Vector256.Create(xScalar3);
                    Vector256<double> aVec03 = Vector256.LoadUnsafe(ref aHead);
                    Vector256<double> aVec13 = Vector256.LoadUnsafe(ref aHead1);
                    xHead = ref Unsafe.Add(ref xHead, xStride);
                    aHead = ref Unsafe.Add(ref aHead, aRowStride);
                    aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                    yVec0 = Fma.MultiplyAdd(aVec03, xVec3, yVec0);
                    yVec1 = Fma.MultiplyAdd(aVec13, xVec3, yVec1);
                }
                for (; i < rows; i++)
                {
                    double xScalar = xHead * alpha;
                    Vector256<double> xVec = Vector256.Create(xScalar);
                    Vector256<double> aVec0 = Vector256.LoadUnsafe(ref aHead);
                    Vector256<double> aVec1 = Vector256.LoadUnsafe(ref aHead1);
                    xHead = ref Unsafe.Add(ref xHead, xStride);
                    aHead = ref Unsafe.Add(ref aHead, aRowStride);
                    aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                    yVec0 = Fma.MultiplyAdd(aVec0, xVec, yVec0);
                    yVec1 = Fma.MultiplyAdd(aVec1, xVec, yVec1);
                }

                yVec0.StoreUnsafe(ref yBufferHead);
                yVec1.StoreUnsafe(ref yBufferHead, 4);
                for (i = 0; i < DotxF_Kernel_RowMajor_Vector256_PerferredCount; i++)
                {
                    yHead *= beta;
                    yHead += yBufferHead;
                    yHead = ref Unsafe.Add(ref yHead, yStride);
                    yBufferHead = ref Unsafe.Add(ref yBufferHead, 1);
                }
            }

            public const nint DotxF_Kernel_ColMajor_Vector256_PerferredCount = 6;

            public static void DotxF_Kernel_ColMajor_Vector256_Perferred_4p6(nint rows,  double alpha, ref double aHead, nint aColStride, ref double xHead, nint xStride,  double beta, ref double yHead, nint yStride)
            {
                using var xBuffer = new BufferDVectorSpan(ref xHead, rows, xStride, alpha);
                xHead = ref xBuffer.bufferHead;
                xStride = 1;
                ApplyWithV_Kernel<MultiplyOperator<double>, double>(ref yHead, 
                    beta, new(DotxF_Kernel_ColMajor_Vector256_PerferredCount, yStride));
                Span<double> vals = stackalloc double[(int)DotxF_Kernel_ColMajor_Vector256_PerferredCount];
                ref double aRef0 = ref aHead;
                ref double aRef1 = ref Unsafe.Add(ref aHead, aColStride);
                ref double aRef2 = ref Unsafe.Add(ref aHead, 2 * aColStride);
                ref double aRef3 = ref Unsafe.Add(ref aHead, 3 * aColStride);
                ref double aRef4 = ref Unsafe.Add(ref aHead, 4 * aColStride);
                ref double aRef5 = ref Unsafe.Add(ref aHead, 5 * aColStride);

                nint i = 0;
                nint iterStride = 1;
                nint iterSize = iterStride * Vector256<double>.Count;
                if (rows >= iterSize)
                {

                    Vector256<double> yVec0 = Vector256<double>.Zero;
                    Vector256<double> yVec1 = Vector256<double>.Zero;
                    Vector256<double> yVec2 = Vector256<double>.Zero;
                    Vector256<double> yVec3 = Vector256<double>.Zero;
                    Vector256<double> yVec4 = Vector256<double>.Zero;
                    Vector256<double> yVec5 = Vector256<double>.Zero;

                    for (; i <= rows - iterSize; i += iterSize)
                    {
                        Vector256<double> aVec0 = Vector256.LoadUnsafe(ref aRef0);
                        Vector256<double> aVec1 = Vector256.LoadUnsafe(ref aRef1);
                        Vector256<double> aVec2 = Vector256.LoadUnsafe(ref aRef2);
                        Vector256<double> aVec3 = Vector256.LoadUnsafe(ref aRef3);
                        Vector256<double> aVec4 = Vector256.LoadUnsafe(ref aRef4);
                        Vector256<double> aVec5 = Vector256.LoadUnsafe(ref aRef5);

                        Vector256<double> xVec = Vector256.LoadUnsafe(ref xHead);

                        yVec0 = Fma.MultiplyAdd(aVec0, xVec, yVec0);
                        yVec1 = Fma.MultiplyAdd(aVec1, xVec, yVec1);
                        yVec2 = Fma.MultiplyAdd(aVec2, xVec, yVec2);
                        yVec3 = Fma.MultiplyAdd(aVec3, xVec, yVec3);
                        yVec4 = Fma.MultiplyAdd(aVec4, xVec, yVec4);
                        yVec5 = Fma.MultiplyAdd(aVec5, xVec, yVec5);

                        aRef0 = ref Unsafe.Add(ref aRef0, iterSize);
                        aRef1 = ref Unsafe.Add(ref aRef1, iterSize);
                        aRef2 = ref Unsafe.Add(ref aRef2, iterSize);
                        aRef3 = ref Unsafe.Add(ref aRef3, iterSize);
                        aRef4 = ref Unsafe.Add(ref aRef4, iterSize);
                        aRef5 = ref Unsafe.Add(ref aRef5, iterSize);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                    }

                    vals[0] = yVec0.Sum();
                    vals[1] = yVec1.Sum();
                    vals[2] = yVec2.Sum();
                    vals[3] = yVec3.Sum();
                    vals[4] = yVec4.Sum();
                    vals[5] = yVec5.Sum();
                }
                for (; i < rows; i++)
                {
                    double val = xHead;
                    vals[0] += aRef0 * val;
                    vals[1] += aRef1 * val;
                    vals[2] += aRef2 * val;
                    vals[3] += aRef3 * val;
                    vals[4] += aRef4 * val;
                    vals[5] += aRef5 * val;

                    aRef0 = ref Unsafe.Add(ref aRef0, 1);
                    aRef1 = ref Unsafe.Add(ref aRef1, 1);
                    aRef2 = ref Unsafe.Add(ref aRef2, 1);
                    aRef3 = ref Unsafe.Add(ref aRef3, 1);
                    aRef4 = ref Unsafe.Add(ref aRef4, 1);
                    aRef5 = ref Unsafe.Add(ref aRef5, 1);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
                ApplyToV_Kernel<AddOperator<double>>
                (ref vals[0], ref yHead, 
                new(DotxF_Kernel_ColMajor_Vector256_PerferredCount, 1,yStride));
            }

            public static void DotxF_Kernel(nint rows, nint cols,  double alpha, ref double aHead, nint aRowStride, nint aColStride, ref double xHead, nint xStride,  double beta, ref double yHead, nint yStride)
            {
                ApplyWithV_Kernel<MultiplyOperator<double>, double>
                    (ref yHead, beta, new(cols, yStride));
                ref double aColHead = ref aHead;
                for (nint i = 0; i < cols; i++)
                {
                    DotV_Kernel(ref aColHead, ref xHead, new(rows, aRowStride, xStride), out var rho);
                    aColHead = ref Unsafe.Add(ref aColHead, aColStride);
                    yHead += alpha * rho;
                    yHead = ref Unsafe.Add(ref yHead, yStride);
                }
            }
        }

        public static void DotxAxpyF( double alpha, Matrix a,  Vector w,  Vector x, 
             double beta,  Vector y,  Vector z)
        {
            var rows = a.Rows;
            var cols = a.Cols;
            if (z.Length != rows || w.Length != rows || x.Length != cols || y.Length != cols)
                throw new ArgumentException("Matrix and vector dimensions do not match.");
            Details.DotxAxpyF_Impl(rows, cols,
                        alpha,
                        ref a[0, 0], a.RowStride, a.ColStride,
                        ref w[0], w.Stride,
                        ref x[0], x.Stride,
                        beta,
                        ref y[0], y.Stride,
                        ref z[0], z.Stride);
        }

        public static partial class Details
        {
            public static void DotxAxpyF_Impl(nint rows, nint cols, 
                 double alpha, 
                ref double aHead, nint aRowStride, nint aColStride, 
                ref double wHead, nint wStride, 
                ref double xHead, nint xStride, 
                 double beta, 
                ref double yHead, nint yStride, 
                ref double zHead, nint zStride)
            {
                if ((rows == 0) || (cols == 0))
                    return;
                else
                {
                    if (cols == 1)
                    {
                        DotV_Impl(ref aHead, ref xHead, new(rows, aRowStride, wStride), out var rho);
                        yHead = beta * yHead + alpha * rho;
                        ApplyToWithV_Impl<AxpyOperator, double>(ref aHead, 
                            alpha * xHead, ref zHead, new(rows, aRowStride, zStride));
                    }
                    else if (aColStride == 1 && cols == DotxAxpyF_Kernel_RowMajor_Vector256_PerferredCount)
                        DotxAxpyF_Kernel_RowMajor_Vector256_Perferred_4p4(rows,
                            alpha,
                            ref aHead, aRowStride,
                            ref wHead, wStride,
                            ref xHead, xStride,
                            beta,
                            ref yHead, yStride,
                            ref zHead, zStride);
                    else if (aRowStride == 1 && cols == DotxAxpyF_Kernel_ColMajor_Vector256_PerferredCount)
                        DotxAxpyF_Kernel_ColMajor_Vector256_Perferred_4p4(rows,
                            alpha,
                            ref aHead, aColStride,
                            ref wHead, wStride,
                            ref xHead, xStride,
                            beta,
                            ref yHead, yStride,
                            ref zHead, zStride);
                    else
                        DotxAxpyF_Kernel(rows, cols,
                            alpha,
                            ref aHead, aRowStride, aColStride,
                            ref wHead, wStride,
                            ref xHead, xStride,
                            beta,
                            ref yHead, yStride,
                            ref zHead, zStride);
                }
            }

            public const nint DotxAxpyF_Kernel_ColMajor_Vector256_PerferredCount = 4;

            public static void DotxAxpyF_Kernel_ColMajor_Vector256_Perferred_4p4(nint rows, 
                double alpha, 
                ref double aHead, nint aColStride, 
                ref double wHead, nint wStride, 
                ref double xHead, nint xStride, 
                double beta, 
                ref double yHead, nint yStride, 
                ref double zHead, nint zStride)
            {
                nint pref = DotxAxpyF_Kernel_ColMajor_Vector256_PerferredCount;
                nint iterStride = 1;
                nint iterSize = iterStride * Vector256<double>.Count;
                nint i = 0;

                double x0 = xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double x1 = xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double x2 = xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double x3 = xHead;

                Span<double> yBuffer = stackalloc double[(int)pref];
                ref var yBufferHead = ref MemoryMarshal.GetReference(yBuffer);

                using var wBuffer = new BufferDVectorSpan(ref wHead, rows, wStride, shouldCopyBack: false);
                wHead = ref wBuffer.bufferHead;
                wStride = 1;

                using var zBuffer = new BufferDVectorSpan(ref zHead, rows, zStride, shouldCopyBack: true);
                zHead = ref zBuffer.bufferHead;
                zStride = 1;

                i = 0;
                ref double aHead1 = ref Unsafe.Add(ref aHead, aColStride);
                ref double aHead2 = ref Unsafe.Add(ref aHead, 2 * aColStride);
                ref double aHead3 = ref Unsafe.Add(ref aHead, 3 * aColStride);
                Vector256<double> alphaVec = Vector256.Create(alpha);

                if (rows >= iterSize)
                {

                    Vector256<double> xVec0 = Vector256.Create(x0);
                    Vector256<double> xVec1 = Vector256.Create(x1);
                    Vector256<double> xVec2 = Vector256.Create(x2);
                    Vector256<double> xVec3 = Vector256.Create(x3);
                    Vector256<double> yVec0 = Vector256<double>.Zero;
                    Vector256<double> yVec1 = Vector256<double>.Zero;
                    Vector256<double> yVec2 = Vector256<double>.Zero;
                    Vector256<double> yVec3 = Vector256<double>.Zero;

                    for (; i <= rows - iterSize; i += iterSize)
                    {
                        Vector256<double> aVec0 = Vector256.LoadUnsafe(ref aHead);
                        aVec0 *= alphaVec;
                        Vector256<double> aVec1 = Vector256.LoadUnsafe(ref aHead1);
                        aVec1 *= alphaVec;
                        Vector256<double> aVec2 = Vector256.LoadUnsafe(ref aHead2);
                        aVec2 *= alphaVec;
                        Vector256<double> aVec3 = Vector256.LoadUnsafe(ref aHead3);
                        aVec3 *= alphaVec;

                        Vector256<double> wVec = Vector256.LoadUnsafe(ref wHead);

                        yVec0 = Fma.MultiplyAdd(aVec0, wVec, yVec0);
                        yVec1 = Fma.MultiplyAdd(aVec1, wVec, yVec1);
                        yVec2 = Fma.MultiplyAdd(aVec2, wVec, yVec2);
                        yVec3 = Fma.MultiplyAdd(aVec3, wVec, yVec3);

                        Vector256<double> zVec = Vector256.LoadUnsafe(ref zHead);
                        zVec = Fma.MultiplyAdd(aVec0, xVec0, zVec);
                        zVec = Fma.MultiplyAdd(aVec1, xVec1, zVec);
                        zVec = Fma.MultiplyAdd(aVec2, xVec2, zVec);
                        zVec = Fma.MultiplyAdd(aVec3, xVec3, zVec);

                        zVec.StoreUnsafe(ref zHead);

                        aHead = ref Unsafe.Add(ref aHead, iterSize);
                        aHead1 = ref Unsafe.Add(ref aHead1, iterSize);
                        aHead2 = ref Unsafe.Add(ref aHead2, iterSize);
                        aHead3 = ref Unsafe.Add(ref aHead3, iterSize);
                        wHead = ref Unsafe.Add(ref wHead, iterSize);
                        zHead = ref Unsafe.Add(ref zHead, iterSize);
                    }

                    yBufferHead = yVec0.Sum();
                    yBufferHead = ref Unsafe.Add(ref yBufferHead, 1);
                    yBufferHead = yVec1.Sum();
                    yBufferHead = ref Unsafe.Add(ref yBufferHead, 1);
                    yBufferHead = yVec2.Sum();
                    yBufferHead = ref Unsafe.Add(ref yBufferHead, 1);
                    yBufferHead = yVec3.Sum();
                    yBufferHead = ref MemoryMarshal.GetReference(yBuffer);
                }
                if(i < rows) 
                {
                    Vector256<double> yVec = Vector256.LoadUnsafe(ref yBufferHead);
                    Vector256<double> xVec = Vector256.Create(x0, x1, x2, x3);
                    for (; i < rows; i++)
                    {
                        Vector256<double> aVec = Vector256.Create(aHead, aHead1, aHead2, aHead3);
                        Vector256<double> wVec = Vector256.Create(wHead);
                        aVec *= alphaVec;
                        yVec = Fma.MultiplyAdd(aVec, wVec, yVec);
                        Vector256<double> zVec = aVec * xVec;
                        zHead += zVec.Sum();

                        aHead = ref Unsafe.Add(ref aHead, 1);
                        aHead1 = ref Unsafe.Add(ref aHead1, 1);
                        aHead2 = ref Unsafe.Add(ref aHead2, 1);
                        aHead3 = ref Unsafe.Add(ref aHead3, 1);
                        wHead = ref Unsafe.Add(ref wHead, 1);
                        zHead = ref Unsafe.Add(ref zHead, 1);
                    }
                    yVec.StoreUnsafe(ref yBufferHead);
                }
                for (i = 0; i < DotxAxpyF_Kernel_ColMajor_Vector256_PerferredCount; i++)
                {
                    yHead = yBufferHead + beta * yHead;
                    yHead = ref Unsafe.Add(ref yHead, yStride);
                    yBufferHead = ref Unsafe.Add(ref yBufferHead, 1);
                }
            }

            public const nint DotxAxpyF_Kernel_RowMajor_Vector256_PerferredCount = 4;

            public static void DotxAxpyF_Kernel_RowMajor_Vector256_Perferred_4p4(nint rows,
                double alpha,
                ref double aHead, nint aRowStride,
                ref double wHead, nint wStride,
                ref double xHead, nint xStride,
                double beta,
                ref double yHead, nint yStride,
                ref double zHead, nint zStride)
            {
                nint pref = DotxAxpyF_Kernel_ColMajor_Vector256_PerferredCount;
                nint iterStride = 1;
                nint iterSize = iterStride * Vector256<double>.Count;
                nint i = 0;

                Span<double> yBuffer = stackalloc double[(int)pref];
                ref var yBufferHead = ref MemoryMarshal.GetReference(yBuffer);
                Vector256<double> yVec = Vector256<double>.Zero;
                Vector256<double> xVec = Vector256.LoadUnsafe(ref xHead);

                double x0 = xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double x1 = xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double x2 = xHead;
                xHead = ref Unsafe.Add(ref xHead, xStride);
                double x3 = xHead;

                using var zBuffer = new BufferDVectorSpan(ref zHead, rows, zStride, shouldCopyBack: true);
                zHead = ref zBuffer.bufferHead;
                zStride = 1;

                i = 0;
                ref double aHead1 = ref Unsafe.Add(ref aHead, aRowStride);
                ref double aHead2 = ref Unsafe.Add(ref aHead, 2 * aRowStride);
                ref double aHead3 = ref Unsafe.Add(ref aHead, 3 * aRowStride);
                Vector256<double> alphaVec = Vector256.Create(alpha);

                if (rows >= iterSize)
                {
                    
                    ref var zHead1 = ref Unsafe.Add(ref zHead, 1);
                    ref var zHead2 = ref Unsafe.Add(ref zHead, 2);
                    ref var zHead3 = ref Unsafe.Add(ref zHead, 3);

                    for (; i <= rows - iterSize; i += iterSize)
                    {
                        Vector256<double> aVec0 = Vector256.LoadUnsafe(ref aHead);
                        aVec0 *= alphaVec;
                        Vector256<double> aVec1 = Vector256.LoadUnsafe(ref aHead1);
                        aVec1 *= alphaVec;
                        Vector256<double> aVec2 = Vector256.LoadUnsafe(ref aHead2);
                        aVec2 *= alphaVec;
                        Vector256<double> aVec3 = Vector256.LoadUnsafe(ref aHead3);
                        aVec3 *= alphaVec;

                        Vector256<double> wVec0 = Vector256.Create(wHead);
                        wHead = ref Unsafe.Add(ref wHead, wStride);
                        Vector256<double> wVec1 = Vector256.Create(wHead);
                        wHead = ref Unsafe.Add(ref wHead, wStride);
                        Vector256<double> wVec2 = Vector256.Create(wHead);
                        wHead = ref Unsafe.Add(ref wHead, wStride);
                        Vector256<double> wVec3 = Vector256.Create(wHead);
                        wHead = ref Unsafe.Add(ref wHead, wStride);

                        yVec = Fma.MultiplyAdd(aVec0, wVec0, yVec);
                        yVec = Fma.MultiplyAdd(aVec1, wVec1, yVec);
                        yVec = Fma.MultiplyAdd(aVec2, wVec2, yVec);
                        yVec = Fma.MultiplyAdd(aVec3, wVec3, yVec);

                        Vector256<double> zVec0 = aVec0 * xVec;
                        Vector256<double> zVec1 = aVec1 * xVec;
                        Vector256<double> zVec2 = aVec2 * xVec;
                        Vector256<double> zVec3 = aVec3 * xVec;
                        zHead += zVec0.Sum();
                        zHead1 += zVec1.Sum();
                        zHead2 += zVec2.Sum();
                        zHead3 += zVec3.Sum();

                        aHead = ref Unsafe.Add(ref aHead, iterSize * aRowStride);
                        aHead1 = ref Unsafe.Add(ref aHead1, iterSize * aRowStride);
                        aHead2 = ref Unsafe.Add(ref aHead2, iterSize * aRowStride);
                        aHead3 = ref Unsafe.Add(ref aHead3, iterSize * aRowStride);
                        zHead = ref Unsafe.Add(ref zHead, iterSize);
                        zHead1 = ref Unsafe.Add(ref zHead1, iterSize);
                        zHead2 = ref Unsafe.Add(ref zHead2, iterSize);
                        zHead3 = ref Unsafe.Add(ref zHead3, iterSize);
                    }
                }
                if (i < rows)
                {
                    for (; i < rows; i++)
                    {
                        Vector256<double> aVec = Vector256.LoadUnsafe(ref aHead);
                        Vector256<double> wVec = Vector256.Create(wHead);
                        aVec *= alphaVec;
                        yVec = Fma.MultiplyAdd(aVec, wVec, yVec);
                        Vector256<double> zVec = aVec * xVec;
                        zHead += zVec.Sum();

                        aHead = ref Unsafe.Add(ref aHead, aRowStride);
                        aHead1 = ref Unsafe.Add(ref aHead1, aRowStride);
                        aHead2 = ref Unsafe.Add(ref aHead2, aRowStride);
                        aHead3 = ref Unsafe.Add(ref aHead3, aRowStride);
                        wHead = ref Unsafe.Add(ref wHead, 1);
                        zHead = ref Unsafe.Add(ref zHead, 1);
                    }
                }
                yVec.StoreUnsafe(ref yBufferHead);
                for (i = 0; i < DotxAxpyF_Kernel_ColMajor_Vector256_PerferredCount; i++)
                {
                    yHead = yBufferHead + beta * yHead;
                    yHead = ref Unsafe.Add(ref yHead, yStride);
                    yBufferHead = ref Unsafe.Add(ref yBufferHead, 1);
                }
            }

            public static void DotxAxpyF_Kernel(nint rows, nint cols,
                 double alpha,
                ref double aHead, nint aRowStride, nint aColStride,
                ref double wHead, nint wStride,
                ref double xHead, nint xStride,
                 double beta,
                ref double yHead, nint yStride,
                ref double zHead, nint zStride)
            {
                for (nint i = 0; i < cols; i++)
                {
                    DotV_Impl(ref aHead, ref wHead, new(rows, aRowStride, wStride), out var rho);
                    ApplyToWithV_Impl<AxpyOperator, double>(ref aHead, 
                        alpha * xHead, ref zHead, new(rows, aRowStride,zStride));
                    aHead = ref Unsafe.Add(ref aHead, aColStride);
                    xHead = ref Unsafe.Add(ref xHead, xStride);
                    yHead *= beta;
                    yHead += alpha * rho;
                    yHead = ref Unsafe.Add(ref yHead, yStride);
                }
            }
        }
    }
}
