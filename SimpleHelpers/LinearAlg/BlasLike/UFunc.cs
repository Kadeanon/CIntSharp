using NetFabric.Numerics.Tensors;
using SimpleHelpers.Indices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Text;
using System.Threading.Tasks;

namespace SimpleHelpers.LinearAlg
{
    public static partial class BlasLike
    {
        public static void ApplySelfV<TAction, TIn>(Vector x, TIn alpha)
            where TAction : struct, IUnaryOperator<TIn, double>
            where TIn : struct
        {
            Details.ApplySelfV_Impl<TAction, TIn>(ref x[0], alpha, x.Indice);
        }

        public static partial class Details
        {
            public static void ApplySelfV_Impl<TAction, TIn>(ref double xHead,
                 TIn alpha, SingleIndice indice)
            where TAction : struct, IUnaryOperator<TIn, double>
            where TIn : struct
            {
                if (indice.Length == 0)
                    return;
                else if (indice.Stride > 1)
                    ApplySelfV_Kernel<TAction, TIn>
                        (ref xHead, alpha, indice);
                else if (System.Numerics.Vector.IsHardwareAccelerated && TAction.IsVectorizable)
                    ApplySelfV_Kernel_Vector<TAction, TIn>
                        (indice.Length, ref xHead, alpha);
                else
                    ApplySelfV_Kernel_Unit<TAction, TIn>
                        (indice.Length, ref xHead, alpha);
            }

            public static void ApplySelfV_Kernel_Vector<TAction, TIn>(nint length,
                ref double xHead, TIn alpha)
                where TAction : struct, IUnaryOperator<TIn, double>
                where TIn : struct
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 3);
                    Vector<TIn> alphaVec = System.Numerics.Vector.Create(alpha);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> xVec1 = System.Numerics.Vector.LoadUnsafe(ref xHead1);
                        Vector<double> xVec2 = System.Numerics.Vector.LoadUnsafe(ref xHead2);
                        Vector<double> xVec3 = System.Numerics.Vector.LoadUnsafe(ref xHead3);
                        xVec0 = TAction.Invoke(in alphaVec);
                        xVec1 = TAction.Invoke(in alphaVec);
                        xVec2 = TAction.Invoke(in alphaVec);
                        xVec3 = TAction.Invoke(in alphaVec);
                        xVec0.StoreUnsafe(ref xHead);
                        xVec1.StoreUnsafe(ref xHead1);
                        xVec2.StoreUnsafe(ref xHead2);
                        xVec3.StoreUnsafe(ref xHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                    }
                    for (; i <= length - Vector<double>.Count; i += Vector<double>.Count)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        xVec0 = TAction.Invoke(in alphaVec);
                        xVec0.StoreUnsafe(ref xHead);
                        xHead = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    }
                }
                for (; i < length; i++)
                {
                    xHead = TAction.Invoke(alpha);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }

            public static void ApplySelfV_Kernel_Unit<TAction, TIn>(nint length,
                ref double xHead, TIn alpha)
                where TAction : struct, IUnaryOperator<TIn, double>
                where TIn : struct
            {
                nint iterSize = 4;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, 1);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        xHead = TAction.Invoke(alpha);
                        xHead1 = TAction.Invoke(alpha);
                        xHead2 = TAction.Invoke(alpha);
                        xHead3 = TAction.Invoke(alpha);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    xHead = TAction.Invoke(alpha);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }

            public static void ApplySelfV_Kernel<TAction, TIn>(ref double xHead,
                TIn alpha, SingleIndice indice)
                where TAction : struct, IUnaryOperator<TIn, double>
                where TIn : struct
            {
                nint iterSize = 4;
                nint iterStride = iterSize * indice.Stride;
                nint i = 0;
                if (indice.Length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, indice.Stride);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, indice.Stride * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, indice.Stride * 3);
                    for (; i <= indice.Length - iterSize; i += iterSize)
                    {
                        xHead = TAction.Invoke(alpha);
                        xHead1 = TAction.Invoke(alpha);
                        xHead2 = TAction.Invoke(alpha);
                        xHead3 = TAction.Invoke(alpha);
                        xHead = ref Unsafe.Add(ref xHead, iterStride);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterStride);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterStride);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterStride);
                    }
                }
                for (; i < indice.Length; i++)
                {
                    xHead = TAction.Invoke(alpha);
                    xHead = ref Unsafe.Add(ref xHead, indice.Stride);
                }
            }
        }

        public static void ApplyV<TAction>(Vector x)
            where TAction : struct, IUnaryOperator<double, double>
        {
            Details.ApplyV_Impl<TAction>(ref x[0], x.Indice);
        }

        public static partial class Details
        {
            public static void ApplyV_Impl<TAction>(ref double xHead,
                SingleIndice indice)
            where TAction : struct, IUnaryOperator<double, double>
            {
                if (indice.Length == 0)
                    return;
                else if (indice.Stride > 1)
                    ApplyV_Kernel<TAction>(ref xHead, indice);
                else if (System.Numerics.Vector.IsHardwareAccelerated && TAction.IsVectorizable)
                    ApplyV_Kernel_Vector<TAction>(indice.Length, ref xHead);
                else
                    ApplyV_Kernel_Unit<TAction>(indice.Length, ref xHead);
            }

            public static void ApplyV_Kernel_Vector<TAction>(nint length, ref double xHead)
                            where TAction : struct, IUnaryOperator<double, double>
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> xVec1 = System.Numerics.Vector.LoadUnsafe(ref xHead1);
                        Vector<double> xVec2 = System.Numerics.Vector.LoadUnsafe(ref xHead2);
                        Vector<double> xVec3 = System.Numerics.Vector.LoadUnsafe(ref xHead3);
                        xVec0 = TAction.Invoke(in xVec0);
                        xVec1 = TAction.Invoke(in xVec1);
                        xVec2 = TAction.Invoke(in xVec2);
                        xVec3 = TAction.Invoke(in xVec3);
                        xVec0.StoreUnsafe(ref xHead);
                        xVec1.StoreUnsafe(ref xHead1);
                        xVec2.StoreUnsafe(ref xHead2);
                        xVec3.StoreUnsafe(ref xHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                    }
                    for (; i <= length - Vector<double>.Count; i += Vector<double>.Count)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        xVec0 = TAction.Invoke(in xVec0);
                        xVec0.StoreUnsafe(ref xHead);
                        xHead = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    }
                }
                for (; i < length; i++)
                {
                    xHead = TAction.Invoke(xHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }

            public static void ApplyV_Kernel_Unit<TAction>(nint length, ref double xHead)
                            where TAction : struct, IUnaryOperator<double, double>
            {
                nint iterSize = 4;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, 1);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        xHead = TAction.Invoke(xHead);
                        xHead1 = TAction.Invoke(xHead1);
                        xHead2 = TAction.Invoke(xHead2);
                        xHead3 = TAction.Invoke(xHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    xHead = TAction.Invoke(xHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }

            public static void ApplyV_Kernel<TAction>(ref double xHead, SingleIndice indice)
                            where TAction : struct, IUnaryOperator<double, double>
            {
                nint iterSize = 4;
                nint iterStride = iterSize * indice.Stride;
                nint i = 0;
                if (indice.Length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, indice.Stride);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, indice.Stride * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, indice.Stride * 3);
                    for (; i <= indice.Length - iterSize; i += iterSize)
                    {
                        xHead = TAction.Invoke(xHead);
                        xHead1 = TAction.Invoke(xHead1);
                        xHead2 = TAction.Invoke(xHead2);
                        xHead3 = TAction.Invoke(xHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterStride);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterStride);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterStride);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterStride);
                    }
                }
                for (; i < indice.Length; i++)
                {
                    xHead = TAction.Invoke(xHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.Stride);
                }
            }
        }

        public static void ApplyWithV<TUnaryInAction, TIn>(Vector x, TIn alpha)
            where TUnaryInAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
        {
            Details.ApplyWithV_Impl<TUnaryInAction, TIn>(ref x[0],
                alpha, x.Indice);
        }

        public static partial class Details
        {
            public static void ApplyWithV_Impl<TAction, TIn>(ref double xHead,
                TIn alpha, SingleIndice indice)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                if (indice.Length == 0)
                    return;
                else if (indice.Stride > 1)
                    ApplyWithV_Kernel<TAction, TIn>
                        (ref xHead, alpha, indice);
                else if (System.Numerics.Vector.IsHardwareAccelerated && TAction.IsVectorizable)
                    ApplyWithV_Kernel_Vector<TAction, TIn>(indice.Length, ref xHead, alpha);
                else
                    ApplyWithV_Kernel_Unit<TAction, TIn>(indice.Length, ref xHead, alpha);
            }

            public static void ApplyWithV_Kernel_Vector<TAction, TIn>(nint length, ref double xHead, TIn alpha)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 3);
                    Vector<TIn> alphaVec = System.Numerics.Vector.Create(alpha);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> xVec1 = System.Numerics.Vector.LoadUnsafe(ref xHead1);
                        Vector<double> xVec2 = System.Numerics.Vector.LoadUnsafe(ref xHead2);
                        Vector<double> xVec3 = System.Numerics.Vector.LoadUnsafe(ref xHead3);
                        xVec0 = TAction.Invoke(in xVec0, in alphaVec);
                        xVec1 = TAction.Invoke(in xVec1, in alphaVec);
                        xVec2 = TAction.Invoke(in xVec2, in alphaVec);
                        xVec3 = TAction.Invoke(in xVec3, in alphaVec);
                        xVec0.StoreUnsafe(ref xHead);
                        xVec1.StoreUnsafe(ref xHead1);
                        xVec2.StoreUnsafe(ref xHead2);
                        xVec3.StoreUnsafe(ref xHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                    }
                    for (; i <= length - Vector<double>.Count; i += Vector<double>.Count)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        xVec0 = TAction.Invoke(in xVec0, in alphaVec);
                        xVec0.StoreUnsafe(ref xHead);
                        xHead = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    }
                }
                for (; i < length; i++)
                {
                    xHead = TAction.Invoke(xHead, alpha);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }

            public static void ApplyWithV_Kernel_Unit<TAction, TIn>(nint length,
                ref double xHead, TIn alpha)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                nint iterSize = 4;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, 1);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        xHead = TAction.Invoke(xHead, alpha);
                        xHead1 = TAction.Invoke(xHead1, alpha);
                        xHead2 = TAction.Invoke(xHead2, alpha);
                        xHead3 = TAction.Invoke(xHead3, alpha);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    xHead = TAction.Invoke(xHead, alpha);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                }
            }

            public static void ApplyWithV_Kernel<TAction, TIn>(ref double xHead,
                TIn alpha, SingleIndice indice)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                nint iterSize = 4;
                nint iterStride = iterSize * indice.Stride;
                nint i = 0;
                if (indice.Length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, indice.Stride);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, indice.Stride * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, indice.Stride * 3);
                    for (; i <= indice.Length - iterSize; i += iterSize)
                    {
                        xHead = TAction.Invoke(xHead, alpha);
                        xHead1 = TAction.Invoke(xHead1, alpha);
                        xHead2 = TAction.Invoke(xHead2, alpha);
                        xHead3 = TAction.Invoke(xHead3, alpha);
                        xHead = ref Unsafe.Add(ref xHead, iterStride);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterStride);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterStride);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterStride);
                    }
                }
                for (; i < indice.Length; i++)
                {
                    xHead = TAction.Invoke(xHead, alpha);
                    xHead = ref Unsafe.Add(ref xHead, indice.Stride);
                }
            }
        }

        public static void ApplyAssignV<TAction>(Vector x, Vector y)
            where TAction : struct, IUnaryOperator<double, double>
        {
            var indice = CheckIndice(x, y);
            Details.ApplyAssignV_Impl<TAction>(ref x.GetHeadRef(),
                ref y.GetHeadRef(), indice);
        }

        public static partial class Details
        {
            public static void ApplyAssignV_Impl<TAction>(ref double xHead,
                ref double yHead, DoubleIndice indice)
            where TAction : struct, IUnaryOperator<double, double>
            {
                if ((indice.Length == 0))
                    return;
                else if ((indice.AStride > 1) || (indice.BStride > 1))
                    ApplyAssignV_Kernel<TAction>(ref xHead, ref yHead, indice);
                else if (System.Numerics.Vector.IsHardwareAccelerated && TAction.IsVectorizable)
                    ApplyAssignV_Kernel_Vector<TAction>(indice.Length, ref xHead, ref yHead);
                else
                    ApplyAssignV_Kernel_Unit<TAction>(indice.Length, ref xHead, ref yHead);
            }

            public static void ApplyAssignV_Kernel_Vector<TAction>(nint length,
                ref double xHead, ref double yHead)
            where TAction : struct, IUnaryOperator<double, double>
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, Vector<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, Vector<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, Vector<double>.Count * 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> xVec1 = System.Numerics.Vector.LoadUnsafe(ref xHead1);
                        Vector<double> xVec2 = System.Numerics.Vector.LoadUnsafe(ref xHead2);
                        Vector<double> xVec3 = System.Numerics.Vector.LoadUnsafe(ref xHead3);
                        Vector<double> yVec0 = TAction.Invoke(in xVec0);
                        Vector<double> yVec1 = TAction.Invoke(in xVec1);
                        Vector<double> yVec2 = TAction.Invoke(in xVec2);
                        Vector<double> yVec3 = TAction.Invoke(in xVec3);
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
                    for (; i <= length - Vector<double>.Count; i += Vector<double>.Count)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> yVec0 = System.Numerics.Vector.LoadUnsafe(ref yHead);
                        yVec0 = TAction.Invoke(in xVec0);
                        yVec0.StoreUnsafe(ref yHead);
                        xHead = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                        yHead = ref Unsafe.Add(ref yHead, Vector<double>.Count);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = TAction.Invoke(xHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void ApplyAssignV_Kernel_Unit<TAction>(nint length,
                ref double xHead, ref double yHead)
            where TAction : struct, IUnaryOperator<double, double>
            {
                nint iterSize = 4;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, 1);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, 1);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        yHead = TAction.Invoke(xHead);
                        yHead1 = TAction.Invoke(xHead1);
                        yHead2 = TAction.Invoke(xHead2);
                        yHead3 = TAction.Invoke(xHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = TAction.Invoke(xHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void ApplyAssignV_Kernel<TAction>(ref double xHead,
                ref double yHead, DoubleIndice indice)
            where TAction : struct, IUnaryOperator<double, double>
            {
                nint i = 0;
                int iterSize = 4;
                nint xIterStride = iterSize * indice.AStride;
                nint yIterStride = iterSize * indice.BStride;
                if (indice.Length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, indice.AStride);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, indice.AStride * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, indice.AStride * 3);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, indice.BStride);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, indice.BStride * 2);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, indice.BStride * 3);
                    for (; i <= indice.Length - iterSize; i += iterSize)
                    {
                        yHead = TAction.Invoke(xHead);
                        yHead1 = TAction.Invoke(xHead1);
                        yHead2 = TAction.Invoke(xHead2);
                        yHead3 = TAction.Invoke(xHead3);
                        xHead = ref Unsafe.Add(ref xHead, xIterStride);
                        xHead1 = ref Unsafe.Add(ref xHead1, xIterStride);
                        xHead2 = ref Unsafe.Add(ref xHead2, xIterStride);
                        xHead3 = ref Unsafe.Add(ref xHead3, xIterStride);
                        yHead = ref Unsafe.Add(ref yHead, yIterStride);
                        yHead1 = ref Unsafe.Add(ref yHead1, yIterStride);
                        yHead2 = ref Unsafe.Add(ref yHead2, yIterStride);
                        yHead3 = ref Unsafe.Add(ref yHead3, yIterStride);
                    }
                }
                for (; i < indice.Length; i++)
                {
                    yHead = TAction.Invoke(xHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
            }
        }

        public static void ApplyWithAssignV<TAction, TIn>(Vector x, TIn alpha, Vector y)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
        {
            var indice = CheckIndice(x, y);
            Details.ApplyWithAssignV_Impl<TAction, TIn>(ref x.GetHeadRef(),
                alpha, ref y.GetHeadRef(), indice);
        }

        public static partial class Details
        {
            public static void ApplyWithAssignV_Impl<TAction, TIn>(ref double xHead,
                TIn alpha, ref double yHead, DoubleIndice indice)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                if (indice.Length == 0)
                    return;
                else if ((indice.AStride > 1) || (indice.BStride > 1))
                    ApplyWithAssignV_Kernel<TAction, TIn>(ref xHead, alpha, ref yHead, indice);
                else if (System.Numerics.Vector.IsHardwareAccelerated && TAction.IsVectorizable)
                    ApplyWithAssignV_Kernel_Vector<TAction, TIn>(indice.Length, ref xHead, alpha, ref yHead);
                else
                    ApplyWithAssignV_Kernel_Unit<TAction, TIn>(indice.Length, ref xHead, alpha, ref yHead);
            }

            public static void ApplyWithAssignV_Kernel_Vector<TAction, TIn>(nint length,
                ref double xHead, TIn alpha, ref double yHead)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, Vector<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, Vector<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, Vector<double>.Count * 3);
                    Vector<TIn> alphaVec = System.Numerics.Vector.Create(alpha);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> xVec1 = System.Numerics.Vector.LoadUnsafe(ref xHead1);
                        Vector<double> xVec2 = System.Numerics.Vector.LoadUnsafe(ref xHead2);
                        Vector<double> xVec3 = System.Numerics.Vector.LoadUnsafe(ref xHead3);
                        Vector<double> yVec0 = TAction.Invoke(in xVec0, in alphaVec);
                        Vector<double> yVec1 = TAction.Invoke(in xVec1, in alphaVec);
                        Vector<double> yVec2 = TAction.Invoke(in xVec2, in alphaVec);
                        Vector<double> yVec3 = TAction.Invoke(in xVec3, in alphaVec);
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
                    for (; i <= length - Vector<double>.Count; i += Vector<double>.Count)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> yVec0 = System.Numerics.Vector.LoadUnsafe(ref yHead);
                        yVec0 = TAction.Invoke(in xVec0, in alphaVec);
                        yVec0.StoreUnsafe(ref yHead);
                        xHead = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                        yHead = ref Unsafe.Add(ref yHead, Vector<double>.Count);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = TAction.Invoke(xHead, alpha);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void ApplyWithAssignV_Kernel_Unit<TAction, TIn>(nint length,
                ref double xHead, TIn alpha, ref double yHead)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                nint iterSize = 4;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, 1);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, 1);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        yHead = TAction.Invoke(xHead, alpha);
                        yHead1 = TAction.Invoke(xHead1, alpha);
                        yHead2 = TAction.Invoke(xHead2, alpha);
                        yHead3 = TAction.Invoke(xHead3, alpha);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = TAction.Invoke(xHead, alpha);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void ApplyWithAssignV_Kernel<TAction, TIn>(ref double xHead,
                TIn alpha, ref double yHead, DoubleIndice indice)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                nint i = 0;
                int iterSize = 4;
                nint iterAStride = iterSize * indice.AStride;
                nint iterBStride = iterSize * indice.BStride;
                if (indice.Length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, indice.AStride);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, indice.AStride * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, indice.AStride * 3);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, indice.BStride);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, indice.BStride * 2);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, indice.BStride * 3);
                    for (; i <= indice.Length - iterSize; i += iterSize)
                    {
                        yHead = TAction.Invoke(xHead, alpha);
                        yHead1 = TAction.Invoke(xHead1, alpha);
                        yHead2 = TAction.Invoke(xHead2, alpha);
                        yHead3 = TAction.Invoke(xHead3, alpha);
                        xHead = ref Unsafe.Add(ref xHead, iterAStride);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterAStride);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterAStride);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterAStride);
                        yHead = ref Unsafe.Add(ref yHead, iterBStride);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterBStride);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterBStride);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterBStride);
                    }
                }
                for (; i < indice.Length; i++)
                {
                    yHead = TAction.Invoke(xHead, alpha);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
            }
        }

        public static void ApplyToV<TAction>(Vector x, Vector y)
            where TAction : struct, IBinaryOperator<double, double, double>
        {
            var indice = CheckIndice(x, y);
            Details.ApplyToV_Impl<TAction>(ref x.GetHeadRef(),
                ref y.GetHeadRef(), indice);
        }

        public static partial class Details
        {
            public static void ApplyToV_Impl<TAction>(ref double xHead,
                 ref double yHead, DoubleIndice indice)
            where TAction : struct, IBinaryOperator<double, double, double>
            {
                if (indice.Length == 0)
                    return;
                else if ((indice.AStride > 1) || (indice.BStride > 1))
                    ApplyToV_Kernel<TAction>(ref xHead, ref yHead, indice);
                else if (System.Numerics.Vector.IsHardwareAccelerated && TAction.IsVectorizable)
                    ApplyToV_Kernel_Vector<TAction>(indice.Length, ref xHead, ref yHead);
                else
                    ApplyToV_Kernel_Unit<TAction>(indice.Length, ref xHead, ref yHead);
            }

            public static void ApplyToV_Kernel_Vector<TAction>(nint length,
                ref double xHead, ref double yHead)
            where TAction : struct, IBinaryOperator<double, double, double>
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, Vector<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, Vector<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, Vector<double>.Count * 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> yVec0 = System.Numerics.Vector.LoadUnsafe(ref yHead);
                        Vector<double> xVec1 = System.Numerics.Vector.LoadUnsafe(ref xHead1);
                        Vector<double> yVec1 = System.Numerics.Vector.LoadUnsafe(ref yHead1);
                        Vector<double> xVec2 = System.Numerics.Vector.LoadUnsafe(ref xHead2);
                        Vector<double> yVec2 = System.Numerics.Vector.LoadUnsafe(ref yHead2);
                        Vector<double> xVec3 = System.Numerics.Vector.LoadUnsafe(ref xHead3);
                        Vector<double> yVec3 = System.Numerics.Vector.LoadUnsafe(ref yHead3);
                        yVec0 = TAction.Invoke(in xVec0, in yVec0);
                        yVec1 = TAction.Invoke(in xVec1, in yVec1);
                        yVec2 = TAction.Invoke(in xVec2, in yVec2);
                        yVec3 = TAction.Invoke(in xVec3, in yVec3);
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
                    for (; i <= length - Vector<double>.Count; i += Vector<double>.Count)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> yVec0 = System.Numerics.Vector.LoadUnsafe(ref yHead);
                        yVec0 = TAction.Invoke(in xVec0, in yVec0);
                        yVec0.StoreUnsafe(ref yHead);
                        xHead = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                        yHead = ref Unsafe.Add(ref yHead, Vector<double>.Count);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = TAction.Invoke(xHead, yHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void ApplyToV_Kernel_Unit<TAction>(nint length,
                ref double xHead, ref double yHead)
            where TAction : struct, IBinaryOperator<double, double, double>
            {
                nint iterSize = 4;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, 1);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, 1);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        yHead = TAction.Invoke(xHead, yHead);
                        yHead1 = TAction.Invoke(xHead1, yHead1);
                        yHead2 = TAction.Invoke(xHead2, yHead2);
                        yHead3 = TAction.Invoke(xHead3, yHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = TAction.Invoke(xHead, yHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void ApplyToV_Kernel<TAction>(ref double xHead,
                ref double yHead, DoubleIndice indice)
            where TAction : struct, IBinaryOperator<double, double, double>
            {
                nint i = 0;
                int iterSize = 4;
                nint iterAStride = iterSize * indice.AStride;
                nint iterBStride = iterSize * indice.BStride;
                if (indice.Length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, indice.AStride);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, indice.AStride * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, indice.AStride * 3);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, indice.BStride);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, indice.BStride * 2);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, indice.BStride * 3);
                    for (; i <= indice.Length - iterSize; i += iterSize)
                    {
                        yHead = TAction.Invoke(xHead, yHead);
                        yHead1 = TAction.Invoke(xHead1, yHead1);
                        yHead2 = TAction.Invoke(xHead2, yHead2);
                        yHead3 = TAction.Invoke(xHead3, yHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterAStride);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterAStride);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterAStride);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterAStride);
                        yHead = ref Unsafe.Add(ref yHead, iterBStride);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterBStride);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterBStride);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterBStride);
                    }
                }
                for (; i < indice.Length; i++)
                {
                    yHead = TAction.Invoke(xHead, yHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
            }
        }

        public static void ApplyToWithV<TAction, TIn>(Vector x, TIn alpha, Vector y)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
        {
            var indice = CheckIndice(x, y);
            Details.ApplyToWithV_Impl<TAction, TIn>(ref x.GetHeadRef(),
                alpha, ref y.GetHeadRef(), indice);
        }

        public static partial class Details
        {
            public static void ApplyToWithV_Impl<TAction, TIn>(ref double xHead,
                TIn alpha, ref double yHead, DoubleIndice indice)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
            {
                if (indice.Length == 0)
                    return;
                else if ((indice.AStride > 1) || (indice.BStride > 1))
                    ApplyToWithV_Kernel<TAction, TIn>(ref xHead, alpha, ref yHead, indice);
                else if (System.Numerics.Vector.IsHardwareAccelerated && TAction.IsVectorizable)
                    ApplyToWithV_Kernel_Vector<TAction, TIn>(indice.Length, ref xHead, alpha, ref yHead);
                else
                    ApplyToWithV_Kernel_Unit<TAction, TIn>(indice.Length, ref xHead, alpha, ref yHead);
            }

            public static void ApplyToWithV_Kernel_Vector<TAction, TIn>(nint length,
                ref double xHead, TIn alpha, ref double yHead)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
            {
                nint iterStride = 4;
                nint iterSize = iterStride * Vector<double>.Count;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, Vector<double>.Count);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, Vector<double>.Count * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, Vector<double>.Count * 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, Vector<double>.Count * 3);
                    Vector<TIn> alphaVec = System.Numerics.Vector.Create(alpha);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> yVec0 = System.Numerics.Vector.LoadUnsafe(ref yHead);
                        Vector<double> xVec1 = System.Numerics.Vector.LoadUnsafe(ref xHead1);
                        Vector<double> yVec1 = System.Numerics.Vector.LoadUnsafe(ref yHead1);
                        Vector<double> xVec2 = System.Numerics.Vector.LoadUnsafe(ref xHead2);
                        Vector<double> yVec2 = System.Numerics.Vector.LoadUnsafe(ref yHead2);
                        Vector<double> xVec3 = System.Numerics.Vector.LoadUnsafe(ref xHead3);
                        Vector<double> yVec3 = System.Numerics.Vector.LoadUnsafe(ref yHead3);
                        yVec0 = TAction.Invoke(in xVec0, in alphaVec, in yVec0);
                        yVec1 = TAction.Invoke(in xVec1, in alphaVec, in yVec1);
                        yVec2 = TAction.Invoke(in xVec2, in alphaVec, in yVec2);
                        yVec3 = TAction.Invoke(in xVec3, in alphaVec, in yVec3);
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
                    for (; i <= length - Vector<double>.Count; i += Vector<double>.Count)
                    {
                        Vector<double> xVec0 = System.Numerics.Vector.LoadUnsafe(ref xHead);
                        Vector<double> yVec0 = System.Numerics.Vector.LoadUnsafe(ref yHead);
                        yVec0 = TAction.Invoke(in xVec0, in alphaVec, in yVec0);
                        yVec0.StoreUnsafe(ref yHead);
                        xHead = ref Unsafe.Add(ref xHead, Vector<double>.Count);
                        yHead = ref Unsafe.Add(ref yHead, Vector<double>.Count);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = TAction.Invoke(xHead, alpha, yHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void ApplyToWithV_Kernel_Unit<TAction, TIn>(nint length,
                ref double xHead, TIn alpha, ref double yHead)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
            {
                nint iterSize = 4;
                nint i = 0;
                if (length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, 1);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, 1);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, 2);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, 3);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, 3);
                    for (; i <= length - iterSize; i += iterSize)
                    {
                        yHead = TAction.Invoke(xHead, alpha, yHead);
                        yHead1 = TAction.Invoke(xHead1, alpha, yHead1);
                        yHead2 = TAction.Invoke(xHead2, alpha, yHead2);
                        yHead3 = TAction.Invoke(xHead3, alpha, yHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterSize);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterSize);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterSize);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterSize);
                        yHead = ref Unsafe.Add(ref yHead, iterSize);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterSize);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterSize);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterSize);
                    }
                }
                for (; i < length; i++)
                {
                    yHead = TAction.Invoke(xHead, alpha, yHead);
                    xHead = ref Unsafe.Add(ref xHead, 1);
                    yHead = ref Unsafe.Add(ref yHead, 1);
                }
            }

            public static void ApplyToWithV_Kernel<TAction, TIn>(ref double xHead,
                TIn alpha, ref double yHead, DoubleIndice indice)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
            {
                nint i = 0;
                int iterSize = 4;
                nint iterAStride = iterSize * indice.AStride;
                nint iterBStride = iterSize * indice.BStride;
                if (indice.Length >= iterSize)
                {
                    ref var xHead1 = ref Unsafe.Add(ref xHead, indice.AStride);
                    ref var xHead2 = ref Unsafe.Add(ref xHead, indice.AStride * 2);
                    ref var xHead3 = ref Unsafe.Add(ref xHead, indice.AStride * 3);
                    ref var yHead1 = ref Unsafe.Add(ref yHead, indice.BStride);
                    ref var yHead2 = ref Unsafe.Add(ref yHead, indice.BStride * 2);
                    ref var yHead3 = ref Unsafe.Add(ref yHead, indice.BStride * 3);
                    for (; i <= indice.Length - iterSize; i += iterSize)
                    {
                        yHead = TAction.Invoke(xHead, alpha, yHead);
                        yHead1 = TAction.Invoke(xHead1, alpha, yHead1);
                        yHead2 = TAction.Invoke(xHead2, alpha, yHead2);
                        yHead3 = TAction.Invoke(xHead3, alpha, yHead3);
                        xHead = ref Unsafe.Add(ref xHead, iterAStride);
                        xHead1 = ref Unsafe.Add(ref xHead1, iterAStride);
                        xHead2 = ref Unsafe.Add(ref xHead2, iterAStride);
                        xHead3 = ref Unsafe.Add(ref xHead3, iterAStride);
                        yHead = ref Unsafe.Add(ref yHead, iterBStride);
                        yHead1 = ref Unsafe.Add(ref yHead1, iterBStride);
                        yHead2 = ref Unsafe.Add(ref yHead2, iterBStride);
                        yHead3 = ref Unsafe.Add(ref yHead3, iterBStride);
                    }
                }
                for (; i < indice.Length; i++)
                {
                    yHead = TAction.Invoke(xHead, alpha, yHead);
                    xHead = ref Unsafe.Add(ref xHead, indice.AStride);
                    yHead = ref Unsafe.Add(ref yHead, indice.BStride);
                }
            }
        }

        public static void ApplySelfM<TAction, TIn>(Matrix a, TIn alpha)
            where TAction : struct, IUnaryOperator<TIn, double>
            where TIn : struct
        {
            var rowIndice = a.RowIndice;
            var colIndice = a.ColIndice;
            if ((rowIndice.Length == 0) || (colIndice.Length == 0))
                return;
            Details.ApplySelfM_Impl<TAction, TIn>(ref a.GetHeadRef(),
                alpha, a.RowIndice, a.ColIndice);
        }

        public static partial class Details
        {
            public static void ApplySelfM_Impl<TAction, TIn>(ref double aHead,
                TIn alpha, SingleIndice rowIndice, SingleIndice colIndice)
            where TAction : struct, IUnaryOperator<TIn, double>
            where TIn : struct
            {
                if (rowIndice.Stride < colIndice.Stride)
                {
                    (rowIndice, colIndice) =
                        (colIndice, rowIndice);
                }
                if (rowIndice.Length == 1)
                    Details.ApplySelfV_Impl<TAction, TIn>(ref aHead, alpha, colIndice);
                else if (colIndice.Length == 1)
                    Details.ApplySelfV_Impl<TAction, TIn>(ref aHead, alpha, rowIndice);
                else if (colIndice.Stride > 1)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplySelfV_Kernel<TAction, TIn>(ref aHead, alpha, colIndice);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.Stride);
                    }
                }
                else if (System.Numerics.Vector.IsHardwareAccelerated &&
                    TAction.IsVectorizable)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplySelfV_Kernel_Vector<TAction, TIn>(
                            colIndice.Length, ref aHead, alpha);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.Stride);
                    }
                }
                else
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplySelfV_Kernel_Unit<TAction, TIn>(colIndice.Length, 
                            ref aHead, alpha);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.Stride);
                    }
                }
            }
        }

        public static void ApplyM<TAction>(Matrix a)
            where TAction : struct, IUnaryOperator<double, double>
        {
            var rowIndice = a.RowIndice;
            var colIndice = a.ColIndice;
            if ((rowIndice.Length == 0) || (colIndice.Length == 0))
                return;
            Details.ApplyM_Impl<TAction>(ref a.GetHeadRef(),
                a.RowIndice, a.ColIndice);
        }

        public static partial class Details
        {
            public static void ApplyM_Impl<TAction>(ref double aHead,
                SingleIndice rowIndice, SingleIndice colIndice)
            where TAction : struct, IUnaryOperator<double, double>
            {
                if (rowIndice.Stride < colIndice.Stride)
                {
                    (rowIndice, colIndice) =
                        (colIndice, rowIndice);
                }
                if (rowIndice.Length == 1)
                    Details.ApplyV_Impl<TAction>(ref aHead, colIndice);
                else if (colIndice.Length == 1)
                    Details.ApplyV_Impl<TAction>(ref aHead, rowIndice);
                else if (colIndice.Stride > 1)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyV_Kernel<TAction>(ref aHead, colIndice);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.Stride);
                    }
                }
                else if (System.Numerics.Vector.IsHardwareAccelerated &&
                    TAction.IsVectorizable)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyV_Kernel_Vector<TAction>(colIndice.Length, ref aHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.Stride);
                    }
                }
                else
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyV_Kernel_Unit<TAction>(colIndice.Length, ref aHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.Stride);
                    }
                }
            }
        }

        public static void ApplyWithM<TAction, TIn>(Matrix a, TIn alpha)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
        {
            var rowIndice = a.RowIndice;
            var colIndice = a.ColIndice;
            if ((rowIndice.Length == 0) || (colIndice.Length == 0))
                return;
            Details.ApplyWithM_Impl<TAction, TIn>(ref a.GetHeadRef(),
                alpha, a.RowIndice, a.ColIndice);
        }

        public static partial class Details
        {
            public static void ApplyWithM_Impl<TAction, TIn>(ref double xHead,
                TIn alpha, SingleIndice rowIndice, SingleIndice colIndice)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                if (rowIndice.Stride < colIndice.Stride)
                {
                    (rowIndice, colIndice) =
                        (colIndice, rowIndice);
                }
                if (rowIndice.Length == 1)
                    Details.ApplyWithV_Impl<TAction, TIn>(ref xHead, alpha, colIndice);
                else if (colIndice.Length == 1)
                    Details.ApplyWithV_Impl<TAction, TIn>(ref xHead, alpha, rowIndice);
                else if (colIndice.Stride > 1)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyWithV_Kernel<TAction, TIn>(ref xHead, alpha, colIndice);
                        xHead = ref Unsafe.Add(ref xHead, rowIndice.Stride);
                    }
                }
                else if (System.Numerics.Vector.IsHardwareAccelerated &&
                    TAction.IsVectorizable)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyWithV_Kernel_Vector<TAction, TIn>(colIndice.Length, ref xHead, alpha);
                        xHead = ref Unsafe.Add(ref xHead, rowIndice.Stride);
                    }
                }
                else
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyWithV_Kernel_Unit<TAction, TIn>(colIndice.Length, ref xHead, alpha);
                        xHead = ref Unsafe.Add(ref xHead, rowIndice.Stride);
                    }
                }
            }
        }

        public static void ApplyAssignM<TAction>(Matrix a, Matrix b)
            where TAction : struct, IUnaryOperator<double, double>
        {
            var (rowIndice, colIndice) = CheckIndice(a, b);
            if ((rowIndice.Length == 0) || (colIndice.Length == 0))
                return;
            Details.ApplyAssignM_Impl<TAction>(ref a.GetHeadRef(), 
                ref b.GetHeadRef(), rowIndice, colIndice);
        }

        public static partial class Details
        {
            public static void ApplyAssignM_Impl<TAction>(ref double aHead, 
            ref double bHead, DoubleIndice rowIndice, DoubleIndice colIndice)
            where TAction : struct, IUnaryOperator<double, double>
            {
                if (rowIndice.BStride < colIndice.BStride)
                {
                    (rowIndice, colIndice) =
                        (colIndice, rowIndice);
                }
                if (rowIndice.Length == 1)
                    Details.ApplyAssignV_Impl<TAction>(ref aHead, ref bHead, colIndice);
                else if (colIndice.Length == 1)
                    Details.ApplyAssignV_Impl<TAction>(ref aHead, ref bHead, rowIndice);
                else if (colIndice.AStride > 1 || colIndice.BStride > 1)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyAssignV_Kernel<TAction>(ref aHead, ref bHead, 
                            colIndice);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
                else if (System.Numerics.Vector.IsHardwareAccelerated &&
                    TAction.IsVectorizable)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyAssignV_Kernel_Vector<TAction>(colIndice.Length, ref aHead,
                            ref bHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
                else
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyAssignV_Kernel_Unit<TAction>(colIndice.Length, ref aHead,
                            ref bHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
            }
        }

        public static void ApplyWithAssignM<TAction, TIn>(Matrix a, TIn alpha, Matrix b)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
        {
            var (rowIndice, colIndice) = CheckIndice(a, b);
            if ((rowIndice.Length == 0) || (colIndice.Length == 0))
                return;
            Details.ApplyWithAssignM_Impl<TAction, TIn>(ref a.GetHeadRef(),
                alpha, ref b.GetHeadRef(), rowIndice, colIndice);
        }

        public static partial class Details
        {
            public static void ApplyWithAssignM_Impl<TAction, TIn>(ref double aHead,
            TIn alpha, ref double bHead, DoubleIndice rowIndice, DoubleIndice colIndice)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
            {
                if (rowIndice.BStride < colIndice.BStride)
                {
                    (rowIndice, colIndice) =
                        (colIndice, rowIndice);
                }
                if (rowIndice.Length == 1)
                    Details.ApplyWithAssignV_Impl<TAction, TIn>(ref aHead,
                        alpha, ref bHead, colIndice);
                else if (colIndice.Length == 1)
                    Details.ApplyWithAssignV_Impl<TAction, TIn>(ref aHead,
                        alpha, ref bHead, rowIndice);
                else if (colIndice.AStride > 1 || colIndice.BStride > 1)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyWithAssignV_Kernel<TAction, TIn>(ref aHead,
                            alpha, ref bHead, colIndice);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
                else if (System.Numerics.Vector.IsHardwareAccelerated &&
                    TAction.IsVectorizable)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyWithAssignV_Kernel_Vector<TAction, TIn>(colIndice.Length,
                            ref aHead, alpha, ref bHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
                else
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyWithAssignV_Kernel_Unit<TAction, TIn>(colIndice.Length, ref aHead, alpha,
                            ref bHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
            }
        }

        public static void ApplyToM<TAction>(Matrix a, Matrix b)
            where TAction : struct, IBinaryOperator<double, double, double>
        {
            var (rowIndice, colIndice) = CheckIndice(a, b);
            if ((rowIndice.Length == 0) || (colIndice.Length == 0))
                return;
            Details.ApplyToM_Impl<TAction>(ref a.GetHeadRef(),
                ref b.GetHeadRef(), rowIndice, colIndice);
        }

        public static partial class Details
        {
            public static void ApplyToM_Impl<TAction>(ref double aHead,
            ref double bHead, DoubleIndice rowIndice, DoubleIndice colIndice)
                where TAction : struct, IBinaryOperator<double, double, double>
            {
                if (rowIndice.BStride < colIndice.BStride)
                {
                    (rowIndice, colIndice) =
                        (colIndice, rowIndice);
                }
                if (rowIndice.Length == 1)
                    Details.ApplyToV_Impl<TAction>(ref aHead, ref bHead, colIndice);
                else if (colIndice.Length == 1)
                    Details.ApplyToV_Impl<TAction>(ref aHead, ref bHead, rowIndice);
                else if (colIndice.AStride > 1 || colIndice.BStride > 1)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyToV_Kernel<TAction>(ref aHead, ref bHead, 
                            colIndice);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
                else if (System.Numerics.Vector.IsHardwareAccelerated &&
                    TAction.IsVectorizable)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyToV_Kernel_Vector<TAction>(colIndice.Length, 
                            ref aHead, ref bHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
                else
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyToV_Kernel_Unit<TAction>(colIndice.Length, 
                            ref aHead, ref bHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
            }

        }

        public static void ApplyToWithM<TAction, TIn>(Matrix a, TIn alpha, Matrix b)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
        {
            var (rowIndice, colIndice) = CheckIndice(a, b);
            if ((rowIndice.Length == 0) || (colIndice.Length == 0))
                return;
            Details.ApplyToWithM_Impl<TAction, TIn>(ref a.GetHeadRef(),
                alpha, ref b.GetHeadRef(), rowIndice, colIndice);
        }

        public static partial class Details
        {
            public static void ApplyToWithM_Impl<TAction, TIn>(ref double aHead,
            TIn alpha, ref double bHead, DoubleIndice rowIndice, DoubleIndice colIndice)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
            {
                if (rowIndice.BStride < colIndice.BStride)
                {
                    (rowIndice, colIndice) =
                        (colIndice, rowIndice);
                }
                if (rowIndice.Length == 1)
                    Details.ApplyToWithV_Impl<TAction, TIn>(ref aHead, 
                        alpha, ref bHead, colIndice);
                else if (colIndice.Length == 1)
                    Details.ApplyToWithV_Impl<TAction, TIn>(ref aHead, 
                        alpha, ref bHead, rowIndice);
                else if (colIndice.AStride > 1 || colIndice.BStride > 1)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyToWithV_Kernel<TAction, TIn>(ref aHead, 
                            alpha, ref bHead, colIndice);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
                else if (System.Numerics.Vector.IsHardwareAccelerated &&
                    TAction.IsVectorizable)
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyToWithV_Kernel_Vector<TAction, TIn>(colIndice.Length,
                            ref aHead, alpha, ref bHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
                else
                {
                    for (nint i = 0; i < rowIndice.Length; i++)
                    {
                        Details.ApplyToWithV_Kernel_Unit<TAction, TIn>(colIndice.Length, ref aHead, alpha,
                            ref bHead);
                        aHead = ref Unsafe.Add(ref aHead, rowIndice.AStride);
                        bHead = ref Unsafe.Add(ref bHead, rowIndice.BStride);
                    }
                }
            }
        }

        private static DoubleIndice CheckIndice(Vector x, Vector y)
        {
            if (x.Length != y.Length)
                throw new ArgumentException("Error: x and y must have the same length.");
            return new(x.Length, x.Stride, y.Stride);
        }

        private static (DoubleIndice rowIndice, DoubleIndice colIndice)
            CheckIndice(Matrix a, Matrix b)
        {
            if (a.Rows != b.Rows || a.Cols != b.Cols)
                throw new ArgumentException("Error: a and b must have the same dimensions.");
            return (new DoubleIndice(a.Rows, a.RowStride, b.RowStride),
                    new DoubleIndice(a.Cols, a.ColStride, b.ColStride));
        }
    }
}
