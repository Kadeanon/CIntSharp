using NetFabric.Numerics.Tensors;
using NetFabric.Numerics.Tensors.Operators;
using SimpleHelpers.Indices;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace SimpleHelpers.LinearAlg
{
    public static partial class BlasLike
    {
        public static void Add(Vector x, Vector y)
            => ApplyToV<AddOperator<double>>(x, y);

        public static void AMax(Vector x, out nint index)
        {
            Details.AMaxV_Kernel(x.Length, ref x[0], x.Stride, out index);
        }

        public static partial class Details
        {
            public static void AMaxV_Impl(nint length, ref double xHead, nint xStride, out nint index)
            {
                if (length == 0)
                    index = 0;
                else // TODO: Use SIMD to accelerate this
                    AMaxV_Kernel(length, ref xHead, xStride, out index);
            }

            [Obsolete("Vectorized AMaxV is not implemented yet.", true)]
            public static void AMaxV_Kernel_Vector256(nint length, ref double xHead, out nint index)
            {
                index = -1;
                throw new NotImplementedException("Vectorized AMaxV is not implemented yet.");
            }

            public static double AMaxV_Kernel(nint length, ref double xHead, nint xStride, out nint index)
            {
                index = -1;
                double max = 0.0;
                for (nint i = 0; i < length; i++)
                {
                    double val = double.Abs(xHead);
                    if ((max < val) || (!double.IsNaN(max) && double.IsNaN(val)))
                    {
                        max = val;
                        index = i;
                    }
                    xHead = ref Unsafe.Add(ref xHead, xStride);
                }
                return max;
            }
        }

        public static void Axpy(double alpha, Vector x, Vector y)
            => ApplyToWithV<Details.AxpyOperator, double>(x, alpha, y);

        public static partial class Details
        {
            public struct AxpyOperator : ITernaryOperator<double, double, double, double>
            {
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                public static double Invoke(double x, double a, double y)
                    => x * a + y;

                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                public static Vector<double> Invoke(ref readonly Vector<double> x,
                    ref readonly Vector<double> a, ref readonly Vector<double> y)
                {
                    if (Fma.IsSupported && Vector<double>.Count == 2)
                    {
                        return Fma.MultiplyAdd(
                            x.AsVector128(), a.AsVector128(), y.AsVector128()).AsVector();
                    }
                    else if (Fma.IsSupported && Vector<double>.Count == 4)
                    {
                        return Fma.MultiplyAdd(
                            x.AsVector256(), a.AsVector256(), y.AsVector256()).AsVector();
                    }
                    else
                    {
                        return x * a + y;
                    }
                }
            }
        }

        public static void Axpby(double alpha, Vector x, double beta, Vector y)
        {
            var indice = CheckIndice(x, y);
            Details.AxpbyV_Impl(ref x[0], alpha, beta, ref y[0], indice);
        }

        public static partial class Details
        {
            public static void AxpbyV_Impl(ref double xHead, double alpha, double beta, 
                ref double yHead, DoubleIndice indice)
            {
                if ((indice.Length == 0) || (alpha == 0.0 && beta == 1.0))
                    return;
                else if (alpha == 0.0)
                    ApplyWithV_Impl<MultiplyOperator<double>, double>
                        (ref yHead, beta, indice.B);
                else if (beta == 1.0)
                    ApplyToWithV_Impl<AxpyOperator, double>
                        (ref xHead, alpha, ref yHead, indice);
                else if ((indice.AStride == 1) && (indice.BStride == 1))
                    AxpbyV_Kernel_Vector256(indice.Length, ref xHead, alpha, beta, ref yHead);
                else
                    AxpbyV_Kernel(ref xHead, alpha, beta, ref yHead, indice);
            }

            public static void AxpbyV_Kernel_Vector256(nint length, ref double xHead, double alpha, double beta, ref double yHead)
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
                    Vector256<double> betaVec = Vector256.Create(beta);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector256<double> xVec0 = Vector256.LoadUnsafe(ref xHead);
                        Vector256<double> yVec0 = Vector256.LoadUnsafe(ref yHead);
                        Vector256<double> xVec1 = Vector256.LoadUnsafe(ref xHead1);
                        Vector256<double> yVec1 = Vector256.LoadUnsafe(ref yHead1);
                        Vector256<double> xVec2 = Vector256.LoadUnsafe(ref xHead2);
                        Vector256<double> yVec2 = Vector256.LoadUnsafe(ref yHead2);
                        Vector256<double> xVec3 = Vector256.LoadUnsafe(ref xHead3);
                        Vector256<double> yVec3 = Vector256.LoadUnsafe(ref yHead3);
                        yVec0 *= betaVec;
                        yVec1 *= betaVec;
                        yVec2 *= betaVec;
                        yVec3 *= betaVec;
                        if (Fma.IsSupported)
                        {
                            yVec0 = Fma.MultiplyAdd(xVec0, alphaVec, yVec0);
                            yVec1 = Fma.MultiplyAdd(xVec1, alphaVec, yVec1);
                            yVec2 = Fma.MultiplyAdd(xVec2, alphaVec, yVec2);
                            yVec3 = Fma.MultiplyAdd(xVec3, alphaVec, yVec3);
                        }
                        else
                        {
                            yVec0 += xVec0 * alphaVec;
                            yVec1 += xVec1 * alphaVec;
                            yVec2 += xVec2 * alphaVec;
                            yVec3 += xVec3 * alphaVec;
                        }
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                        yVec0.StoreUnsafe(ref yHead);
                        yVec1.StoreUnsafe(ref yHead1);
                        yVec2.StoreUnsafe(ref yHead2);
                        yVec3.StoreUnsafe(ref yHead3);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    yHead *= beta;
                    yHead += alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void AxpbyV_Kernel(ref double xHead, double alpha, double beta, ref double yHead, DoubleIndice indice)
            {
                nint i = 0;
                for (; i <= indice.Length - 4; i += 4)
                {
                    yHead *= beta;
                    yHead += alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    yHead *= beta;
                    yHead += alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    yHead *= beta;
                    yHead += alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    yHead *= beta;
                    yHead += alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
                for (; i < indice.Length; i++)
                {
                    yHead *= beta;
                    yHead += alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
            }
        }

        public static void Copy(Vector x, Vector y)
            => ApplyAssignV<IdentityOperator<double>>(x, y);

        public static double Dot(Vector x, Vector y)
        {
            var indice = CheckIndice(x, y);
            Details.DotV_Impl(ref x[0], ref y[0], indice, out var rho);
            return rho;
        }

        public static partial class Details
        {
            public static void DotV_Impl(ref double xHead, ref double yHead, DoubleIndice indice, out double rho)
            {
                if (indice.Length == 0)
                    rho = 0.0;
                else if ((indice.AStride == 1) && (indice.BStride == 1))
                    DotV_Kernel_Vector256(indice.Length, ref xHead, ref yHead, out rho);
                else
                    DotV_Kernel(ref xHead, ref yHead, indice, out rho);
            }

            private static void DotV_Kernel_Vector256(nint length, ref double xHead, ref double yHead, out double rho)
            {
                double rho2 = 0.0;
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

                    Vector256<double> rhoVec0 = Vector256<double>.Zero;
                    Vector256<double> rhoVec1 = Vector256<double>.Zero;
                    Vector256<double> rhoVec2 = Vector256<double>.Zero;
                    Vector256<double> rhoVec3 = Vector256<double>.Zero;
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector256<double> xVec0 = Vector256.LoadUnsafe(ref xHead);
                        Vector256<double> yVec0 = Vector256.LoadUnsafe(ref yHead);
                        Vector256<double> xVec1 = Vector256.LoadUnsafe(ref xHead1);
                        Vector256<double> yVec1 = Vector256.LoadUnsafe(ref yHead1);
                        Vector256<double> xVec2 = Vector256.LoadUnsafe(ref xHead2);
                        Vector256<double> yVec2 = Vector256.LoadUnsafe(ref yHead2);
                        Vector256<double> xVec3 = Vector256.LoadUnsafe(ref xHead3);
                        Vector256<double> yVec3 = Vector256.LoadUnsafe(ref yHead3);
                        if (Fma.IsSupported)
                        {
                            rhoVec0 = Fma.MultiplyAdd(xVec0, yVec0, rhoVec0);
                            rhoVec1 = Fma.MultiplyAdd(xVec1, yVec1, rhoVec1);
                            rhoVec2 = Fma.MultiplyAdd(xVec2, yVec2, rhoVec2);
                            rhoVec3 = Fma.MultiplyAdd(xVec3, yVec3, rhoVec3);
                        }
                        else
                        {
                            rhoVec0 += xVec0 * yVec0;
                            rhoVec1 += xVec1 * yVec1;
                            rhoVec2 += xVec2 * yVec2;
                            rhoVec3 += xVec3 * yVec3;
                        }
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                    }
                    rhoVec0 += rhoVec1;
                    rhoVec2 += rhoVec3;
                    rhoVec0 += rhoVec2;
                    rho2 += rhoVec0.Sum();
                }
                for (; i < length; i++)
                {
                    rho2 += xHead * yHead;
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
                rho = rho2;
            }

            public static void DotV_Kernel(ref double xHead, ref double yHead, DoubleIndice indice, out double rho)
            {
                rho = 0.0; 
                nint i = 0;
                for (; i < indice.Length - 4; i += 4)
                {
                    rho += xHead * yHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    rho += xHead * yHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    rho += xHead * yHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    rho += xHead * yHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
                for (; i < indice.Length; i++)
                {
                    rho += xHead * yHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
            }
        }

        public static void Invert(Vector x)
            => ApplyWithV<Details.InvertOperator, double>(x, 1);

        public static partial class Details
        {
            public struct InvertOperator : IBinaryOperator<double, double, double>
            {
                public static double Invoke(double x, double one)
                    => one / x;

                public static Vector<double> Invoke(
                    ref readonly Vector<double> x, ref readonly Vector<double> ones)
                    => ones / x;

                static bool IOperator.IsVectorizable => true;
            }
        }

        public static void InvScal(double alpha, Vector x)
        {
            if (alpha == 0.0 || alpha == -0.0 || !double.IsFinite(alpha))
            {
                throw new ArgumentException("Error: alpha must be a finite non-zero value.");
            }
            ApplyWithV<MultiplyOperator<double>, double>(x, 1 / alpha);
        }

        public static void Scal(double alpha, Vector x)
            => ApplyWithV<MultiplyOperator<double>, double>(x, alpha);

        public static void Scal2(double alpha, Vector x, Vector y)
            => ApplyWithAssignV<MultiplyOperator<double>, double>(x, alpha, y);

        public static void Pow(double alpha, Vector x)
            => ApplyWithV<Details.PowOperator, double>(x, alpha);

        public static partial class Details
        {
            public struct PowOperator : IBinaryOperator<double, double, double>
            {
                public static double Invoke(double x, double y)
                    => Math.Pow(x, y);

                public static Vector<double> Invoke(ref readonly Vector<double> x, ref readonly Vector<double> y)
                    => throw new NotSupportedException("Vectorized PowV is not supported.");

                static bool IOperator.IsVectorizable => false;
            }
        }

        public static void Sqrt(Vector x)
            => ApplyV<SqrtOperator<double>>(x);

        public static void Set(double alpha, Vector x)
            => ApplySelfV<IdentityOperator<double>, double>(x, alpha);

        public static void Sub(Vector x, Vector y)
            => ApplyToV<Details.SubtractedByOperator>(x, y);

        public static partial class Details
        {
            public struct SubtractedByOperator : IBinaryOperator<double, double, double>
            {
                public static double Invoke(double other, double self)
                    => self - other;

                public static Vector<double> Invoke(
                    ref readonly Vector<double> other, ref readonly Vector<double> self)
                    => self - other;

                static bool IOperator.IsVectorizable => true;
            }
        }

        public static void Swap(Vector x, Vector y)
        {
            var indice = CheckIndice(x, y);
            Details.SwapV_Impl(ref x[0], ref y[0], indice);
        }

        public static partial class Details
        {
            public static void SwapV_Impl(ref double xHead, 
                ref double yHead, DoubleIndice indice)
            {
                if (indice.Length == 0)
                    return;
                else if ((indice.AStride == 1) && (indice.BStride == 1))
                    SwapV_Kernel_Vector256(indice.Length, ref xHead, ref yHead);
                else
                    SwapV_Kernel(ref xHead, ref yHead, indice);
            }

            private static void SwapV_Kernel_Vector256(nint length, ref double xHead, ref double yHead)
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector256<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, 4);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, 4);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, 8);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, 8);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, 12);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, 12);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector256<double> xVec0 = Vector256.LoadUnsafe(ref xHead);
                        Vector256<double> yVec0 = Vector256.LoadUnsafe(ref yHead);
                        Vector256<double> xVec1 = Vector256.LoadUnsafe(ref xHead1);
                        Vector256<double> yVec1 = Vector256.LoadUnsafe(ref yHead1);
                        Vector256<double> xVec2 = Vector256.LoadUnsafe(ref xHead2);
                        Vector256<double> yVec2 = Vector256.LoadUnsafe(ref yHead2);
                        Vector256<double> xVec3 = Vector256.LoadUnsafe(ref xHead3);
                        Vector256<double> yVec3 = Vector256.LoadUnsafe(ref yHead3);
                        Vector256.StoreUnsafe(xVec0, ref yHead);
                        Vector256.StoreUnsafe(xVec1, ref yHead1);
                        Vector256.StoreUnsafe(xVec2, ref yHead2);
                        Vector256.StoreUnsafe(xVec3, ref yHead3);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                        Vector256.StoreUnsafe(yVec0, ref xHead);
                        Vector256.StoreUnsafe(yVec1, ref xHead1);
                        Vector256.StoreUnsafe(yVec2, ref xHead2);
                        Vector256.StoreUnsafe(yVec3, ref xHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    (xHead, yHead) = (yHead, xHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void SwapV_Kernel(ref double xHead, ref double yHead, DoubleIndice indice)
            {
                nint i = 0;
                for (; i < indice.Length - 4; i += 4)
                {
                    (xHead, yHead) = (yHead, xHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    (xHead, yHead) = (yHead, xHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    (xHead, yHead) = (yHead, xHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    (xHead, yHead) = (yHead, xHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
                for (; i < indice.Length; i++)
                {
                    (xHead, yHead) = (yHead, xHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
            }
        }

        public static void Xpby(Vector x, double beta, Vector y)
            => ApplyToWithV<Details.XpbyOperator, double>(x, beta, y);

        public static partial class Details
        {
            public struct XpbyOperator : ITernaryOperator<double, double, double, double>
            {
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                public static double Invoke(double x, double a, double y)
                    => x + a * y;

                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                public static Vector<double> Invoke(ref readonly Vector<double> x,
                    ref readonly Vector<double> a, ref readonly Vector<double> y)
                {
                    if (Fma.IsSupported && Vector<double>.Count == 2)
                    {
                        return Fma.MultiplyAdd(
                            a.AsVector128(), y.AsVector128(), x.AsVector128()).AsVector();
                    }
                    else if (Fma.IsSupported && Vector<double>.Count == 4)
                    {
                        return Fma.MultiplyAdd(
                            a.AsVector256(), y.AsVector256(), x.AsVector256()).AsVector();
                    }
                    else
                    {
                        return x + a * y;
                    }
                }
            }
        }

        public static void Axpy2(double alphax, double alphay, Vector x, Vector y, Vector z)
        {
            TripleIndice indice = CheckIndice(x, y, z);
            Details.Axpy2V_Impl(ref x[0], alphax, alphay, ref y[0], ref z[0], indice);
        }

        public static partial class Details
        {
            public static void Axpy2V_Impl(ref double xHead, double alphax, double alphay, 
                ref double yHead, ref double zHead, TripleIndice indice)
            {
                if (indice.Length == 0)
                    return;
                else if (alphax == 0.0)
                    ApplyToWithV_Impl<AxpyOperator, double>
                        (ref xHead, alphay, ref zHead, indice.BC);
                else if (alphay == 0.0)
                    ApplyToWithV_Impl<AxpyOperator, double>
                        (ref xHead, alphax, ref zHead, indice.AC);
                else if ((indice.AStride == 1) && (indice.BStride == 1) && (indice.CStride == 1))
                    Axpy2V_Kernel_Vector256(indice.Length, alphax, alphay, ref xHead, 
                        ref yHead, ref zHead);
                else
                    Axpy2V_Kernel(ref xHead, alphax, alphay, ref yHead, ref zHead, indice);
            }

            public static void Axpy2V_Kernel_Vector256(nint length, double alphax, double alphay, ref double xHead, ref double yHead, ref double zHead)
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector256<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector256<double>.Count);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, Vector256<double>.Count);
                    ref var zHead1 = ref Unsafe.Add(ref zHead, Vector256<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector256<double>.Count * 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, Vector256<double>.Count * 2);
                    ref var zHead2 = ref Unsafe.Add(ref zHead, Vector256<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector256<double>.Count * 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, Vector256<double>.Count * 3);
                    ref var zHead3 = ref Unsafe.Add(ref zHead, Vector256<double>.Count * 3);
                    Vector256<double> alphaVecx = Vector256.Create(alphax);
                    Vector256<double> alphaVecy = Vector256.Create(alphay);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector256<double> xVec0 = Vector256.LoadUnsafe(ref xHead);
                        Vector256<double> yVec0 = Vector256.LoadUnsafe(ref yHead);
                        Vector256<double> zVec0 = Vector256.LoadUnsafe(ref zHead);
                        Vector256<double> xVec1 = Vector256.LoadUnsafe(ref xHead1);
                        Vector256<double> yVec1 = Vector256.LoadUnsafe(ref yHead1);
                        Vector256<double> zVec1 = Vector256.LoadUnsafe(ref zHead1);
                        Vector256<double> xVec2 = Vector256.LoadUnsafe(ref xHead2);
                        Vector256<double> yVec2 = Vector256.LoadUnsafe(ref yHead2);
                        Vector256<double> zVec2 = Vector256.LoadUnsafe(ref zHead2);
                        Vector256<double> xVec3 = Vector256.LoadUnsafe(ref xHead3);
                        Vector256<double> yVec3 = Vector256.LoadUnsafe(ref yHead3);
                        Vector256<double> zVec3 = Vector256.LoadUnsafe(ref zHead3);
                        if (Fma.IsSupported)
                        {
                            zVec0 = Fma.MultiplyAdd(xVec0, alphaVecx, zVec0);
                            zVec1 = Fma.MultiplyAdd(xVec1, alphaVecx, zVec1);
                            zVec2 = Fma.MultiplyAdd(xVec2, alphaVecx, zVec2);
                            zVec3 = Fma.MultiplyAdd(xVec3, alphaVecx, zVec3);
                            zVec0 = Fma.MultiplyAdd(yVec0, alphaVecy, zVec0);
                            zVec1 = Fma.MultiplyAdd(yVec1, alphaVecy, zVec1);
                            zVec2 = Fma.MultiplyAdd(yVec2, alphaVecy, zVec2);
                            zVec3 = Fma.MultiplyAdd(yVec3, alphaVecy, zVec3);
                        }
                        else
                        {
                            zVec0 += xVec0 * alphaVecx + yVec0 * alphaVecy;
                            zVec1 += xVec1 * alphaVecx + yVec1 * alphaVecy;
                            zVec2 += xVec2 * alphaVecx + yVec2 * alphaVecy;
                            zVec3 += xVec3 * alphaVecx + yVec3 * alphaVecy;
                        }
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                        zVec0.StoreUnsafe(ref zHead);
                        zVec1.StoreUnsafe(ref zHead1);
                        zVec2.StoreUnsafe(ref zHead2);
                        zVec3.StoreUnsafe(ref zHead3);
                        zHead = ref Unsafe.Add(ref zHead, iterSize);
                        zHead1 = ref Unsafe.Add(ref zHead1, iterSize);
                        zHead2 = ref Unsafe.Add(ref zHead2, iterSize);
                        zHead3 = ref Unsafe.Add(ref zHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    zHead += alphax * xHead;
                    zHead += alphay * yHead;
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                    zHead = ref Unsafe.Add(ref zHead, 1);
                }
            }

            public static void Axpy2V_Kernel(ref double xHead, double alphax, double alphay, ref double yHead, ref double zHead, TripleIndice indice)
            {
                for (nint i = 0; i < indice.Length; i++)
                {
                    double sum = alphax * xHead + alphay * yHead;
                    zHead += sum;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    zHead = ref Unsafe.Add(ref zHead, indice.CStride);
                }
            }
        }

        public static void DotAxpy(double alpha,  Vector x, Vector y, out double rho, Vector z)
        {
            TripleIndice indice = CheckIndice(x, y, z);
            Details.DotAxpyV_Impl(ref x[0], alpha, ref y[0], out rho, ref z[0], indice);
        }

        public static partial class Details
        {
            public static void DotAxpyV_Impl(ref double xHead, double alpha, ref double yHead, 
                out double rho, ref double zHead, TripleIndice indice)
            {
                if (indice.Length == 0)
                    rho = 0.0;
                else if (alpha == 0.0)
                    DotV_Impl(ref xHead, ref yHead, indice.AB, out rho);
                else if ((indice.AStride == 1) && (indice.BStride == 1) && (indice.CStride == 1))
                    DotAxpyV_Kernel_Vector256(indice.Length, alpha, 
                        ref xHead, ref yHead, out rho, ref zHead);
                else
                    DotAxpyV_Kernel(ref xHead, alpha, ref yHead, out rho, ref zHead, indice);
            }

            private static void DotAxpyV_Kernel_Vector256(nint length, double alpha, ref double xHead, ref double yHead, out double rho, ref double zHead)
            {
                rho = 0.0;
                nint iterStride = 3;
                nint iterSize = iterStride * Vector256<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector256<double>.Count);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, Vector256<double>.Count);
                    ref var zHead1 = ref Unsafe.Add(ref zHead, Vector256<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector256<double>.Count * 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, Vector256<double>.Count * 2);
                    ref var zHead2 = ref Unsafe.Add(ref zHead, Vector256<double>.Count * 2);
                    Vector256<double> rhoVec0 = Vector256<double>.Zero;
                    Vector256<double> rhoVec1 = Vector256<double>.Zero;
                    Vector256<double> rhoVec2 = Vector256<double>.Zero;
                    Vector256<double> alphaVec = Vector256.Create(alpha);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector256<double> xVec0 = Vector256.LoadUnsafe(ref xHead);
                        Vector256<double> yVec0 = Vector256.LoadUnsafe(ref yHead);
                        Vector256<double> zVec0 = Vector256.LoadUnsafe(ref zHead);
                        Vector256<double> xVec1 = Vector256.LoadUnsafe(ref xHead1);
                        Vector256<double> yVec1 = Vector256.LoadUnsafe(ref yHead1);
                        Vector256<double> zVec1 = Vector256.LoadUnsafe(ref zHead1);
                        Vector256<double> xVec2 = Vector256.LoadUnsafe(ref xHead2);
                        Vector256<double> yVec2 = Vector256.LoadUnsafe(ref yHead2);
                        Vector256<double> zVec2 = Vector256.LoadUnsafe(ref zHead2);
                        if (Fma.IsSupported)
                        {
                            rhoVec0 = Fma.MultiplyAdd(xVec0, yVec0, rhoVec0);
                            rhoVec1 = Fma.MultiplyAdd(xVec1, yVec1, rhoVec1);
                            rhoVec2 = Fma.MultiplyAdd(xVec2, yVec2, rhoVec2);
                            zVec0 = Fma.MultiplyAdd(xVec0, alphaVec, zVec0);
                            zVec1 = Fma.MultiplyAdd(xVec1, alphaVec, zVec1);
                            zVec2 = Fma.MultiplyAdd(xVec2, alphaVec, zVec2);
                        }
                        else
                        {
                            rhoVec0 += xVec0 * yVec0;
                            rhoVec1 += xVec1 * yVec1;
                            rhoVec2 += xVec2 * yVec2;
                            zVec0 += xVec0 * alphaVec;
                            zVec1 += xVec1 * alphaVec;
                            zVec2 += xVec2 * alphaVec;
                        }
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        zVec0.StoreUnsafe(ref zHead);
                        zVec1.StoreUnsafe(ref zHead1);
                        zVec2.StoreUnsafe(ref zHead2);
                        zHead = ref Unsafe.Add(ref zHead, iterSize);
                        zHead1 = ref Unsafe.Add(ref zHead1, iterSize);
                        zHead2 = ref Unsafe.Add(ref zHead2, iterSize);
                    }
                    rhoVec0 += rhoVec1;
                    rhoVec0 += rhoVec2;
                    rho += rhoVec0.Sum();
                }
                for (; i < length; i++)
                {
                    rho += xHead * yHead;
                    zHead += alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                    zHead = ref Unsafe.Add(ref zHead, 1);
                }
            }

            public static void DotAxpyV_Kernel(ref double xHead, double alpha,
                ref double yHead, out double rho, ref double zHead, TripleIndice indice)
            {
                rho = 0.0;
                for (nint i = 0; i < indice.Length; i++)
                {
                    rho += xHead * yHead;
                    zHead += alpha * xHead;
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                    zHead = ref Unsafe.Add(ref zHead, indice.CStride);
                }
            }
        }

        private static TripleIndice CheckIndice(Vector x, Vector y, Vector z)
        {
            if (x.Length != y.Length)
                throw new ArgumentException("Error: x and y must have the same length.");
            if (x.Length != z.Length)
                throw new ArgumentException("Error: x and z must have the same length.");
            return new(x.Length, x.Stride, y.Stride, z.Stride);
        }

    }
}
