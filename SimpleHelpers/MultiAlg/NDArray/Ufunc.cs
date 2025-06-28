using NetFabric.Numerics.Tensors;
using SimpleHelpers.Indices;
using SimpleHelpers.LinearAlg;
using SimpleHelpers.Utilities;
using System.Runtime.CompilerServices;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {
        public static void ApplySelf<TAction, TIn>(NDArray array, TIn alpha)
            where TAction : struct, IUnaryOperator<TIn, double>
            where TIn : struct
        {
            if (array.Rank == 0)
                return;
            else if(array.Rank == 1)
            {
                BlasLike.Details.ApplySelfV_Impl<TAction, TIn>(ref array.GetHeadRef(),
                    alpha, new(array.Lengths[0], array.Strides[0]));
            }
            else
            {
                Parallel.For(0, array.Lengths[0], index =>
                {
                    ref var head = ref Unsafe.Add(ref array.GetHeadRef(),
                        (nint)index * array.Lengths[0]);
                    Span<SingleIndice> indices = 
                    stackalloc SingleIndice[array.Rank - 1];
                    for (int i = 0; i < array.Rank - 1; i++)
                    {
                        indices[i] = new SingleIndice(array.Lengths[i + 1],
                            array.Strides[i + 1]);
                    }
                    ApplySelf_Impl<TAction, TIn>(ref head, alpha, indices);
                });
            }
        }

        private static void ApplySelf_Impl<TAction, TIn>(ref double head,
            TIn alpha, ReadOnlySpan<SingleIndice> indices)
            where TAction : struct, IUnaryOperator<TIn, double>
            where TIn : struct
        {
            if (indices.Length == 0)
                return;
            else if (indices.Length == 1)
            {
                BlasLike.Details.ApplySelfV_Impl<TAction, TIn>(ref head,
                    alpha, indices[0]);
            }
            else if (indices.Length == 2)
            {
                BlasLike.Details.ApplySelfM_Impl<TAction, TIn>(ref head,
                    alpha, indices[0], indices[1]);
            }
            else
            {
                var firstIndice = indices[0];
                var subIndices = indices[1..];
                for (var i = 0; i < firstIndice.Length; i++)
                {
                    ApplySelf_Impl<TAction, TIn>(ref head, alpha, subIndices);
                    head = ref Unsafe.Add(ref head, firstIndice.Stride);
                }
            }
        }

        public static void BatchApplySelf<TAction>(NDArray array, Vector alpha)
            where TAction : struct, IUnaryOperator<double, double>
        {
            if (array.Rank == 0)
                throw new ArgumentException($"The {nameof(array)} cannot be a scalar array.");
            if (array.Lengths[0] != alpha.Length)
                throw new ArgumentException($"The first dimension of {nameof(array)} " +
                    $"must be equal to the length of {nameof(alpha)}.");

            if (array.Rank == 1)
            {
                BlasLike.ApplyAssignV<TAction>(alpha, array.AsVector());
            }
            else
            {
                Parallel.For(0, array.Lengths[0], index =>
                {
                    var subArray = array.SliceFirstDim((nint)index);
                    ref var head = ref subArray.Data[subArray.Offset];
                    Span<SingleIndice> indices =
                    stackalloc SingleIndice[subArray.Rank];
                    for (int i = 0; i < subArray.Rank; i++)
                    {
                        indices[i] = new SingleIndice(subArray.Lengths[i],
                            subArray.Strides[i]);
                    }
                    ApplySelf_Impl<TAction, double>(ref head, alpha[(nint)index],
                        indices);
                });
            }
        }

        public static void Apply<TAction>(NDArray array)
            where TAction : struct, IUnaryOperator<double, double>
        {
            if (array.Rank == 0)
                return;
            else if (array.Rank == 1)
            {
                BlasLike.Details.ApplyV_Impl<TAction>(ref array.GetHeadRef(),
                    new(array.Lengths[0], array.Strides[0]));
            }
            else
            {
                Parallel.For(0, array.Lengths[0], index =>
                {
                    var subArray = array.SliceFirstDim((nint)index);
                    ref var head = ref subArray.Data[subArray.Offset];
                    Span<SingleIndice> indices =
                    stackalloc SingleIndice[subArray.Rank];
                    for (int i = 0; i < subArray.Rank; i++)
                    {
                        indices[i] = new SingleIndice(subArray.Lengths[i],
                            subArray.Strides[i]);
                    }
                    Apply_Impl<TAction>(ref head, indices);
                });
            }
        }

        private static void Apply_Impl<TAction>(ref double head,
            ReadOnlySpan<SingleIndice> indices)
            where TAction : struct, IUnaryOperator<double, double>
        {
            if (indices.Length == 0)
                return;
            else if (indices.Length == 1)
            {
                BlasLike.Details.ApplyV_Impl<TAction>(ref head, indices[0]);
            }
            else if (indices.Length == 2)
            {
                BlasLike.Details.ApplyM_Impl<TAction>(ref head, 
                    indices[0], indices[1]);
            }
            else
            {
                var firstIndice = indices[0];
                var subIndices = indices[1..];
                for (var i = 0; i < firstIndice.Length; i++)
                {
                    Apply_Impl<TAction>(ref head, subIndices);
                    head = ref Unsafe.Add(ref head, firstIndice.Stride);
                }
            }
        }

        public static void ApplyWith<TAction, TIn>(NDArray array, TIn alpha)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
        {
            if (array.Rank == 0)
                return;
            else if (array.Rank == 1)
            {
                BlasLike.Details.ApplyWithV_Impl<TAction, TIn>(ref array.GetHeadRef(),
                    alpha, new(array.Lengths[0], array.Strides[0]));
            }
            else
            {
                Parallel.For(0, array.Lengths[0], index =>
                {
                    var subArray = array.SliceFirstDim((nint)index);
                    ref var head = ref subArray.Data[subArray.Offset];
                    Span<SingleIndice> indices =
                    stackalloc SingleIndice[subArray.Rank];
                    for (int i = 0; i < subArray.Rank; i++)
                    {
                        indices[i] = new SingleIndice(subArray.Lengths[i],
                            subArray.Strides[i]);
                    }
                    ApplyWith_Impl<TAction, TIn>(ref head, alpha, indices);
                });
            }
        }

        private static void ApplyWith_Impl<TAction, TIn>(ref double head,
            TIn alpha, ReadOnlySpan<SingleIndice> indices)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
        {

            if (indices.Length == 0)
                return;
            else if (indices.Length == 1)
            {
                BlasLike.Details.ApplyWithV_Impl<TAction, TIn>(
                    ref head, alpha, indices[0]);
            }
            else if (indices.Length == 2)
            {
                BlasLike.Details.ApplyWithM_Impl<TAction, TIn>(ref head,
                    alpha, indices[0], indices[1]);
            }
            else
            {
                var firstIndice = indices[0];
                var subIndices = indices[1..];
                for (var i = 0; i < firstIndice.Length; i++)
                {
                    ApplyWith_Impl<TAction, TIn>(ref head, alpha, subIndices);
                    head = ref Unsafe.Add(ref head, firstIndice.Stride);
                }
            }
        }

        public static void BatchApplyWith<TAction>(NDArray array, Vector alpha)
            where TAction : struct, IBinaryOperator<double, double, double>
        {
            if (array.Rank == 0)
                throw new ArgumentException($"The {nameof(array)} cannot be a scalar array.");
            var batchLength = array.Lengths[0];
            if (alpha.Length != batchLength)
                throw new ArgumentException($"The first dimension of {nameof(array)} " +
                    $"must be equal to the length of {nameof(alpha)}.");

            if (array.Rank == 1)
            {
                ref var head = ref array.GetHeadRef();
                ref var alphaHead = ref alpha.GetHeadRef();
                for (nint index = 0; index < array.Lengths[0]; index++)
                {
                    head = TAction.Invoke(head, alphaHead);
                    head = ref Unsafe.Add(ref head, array.Strides[0]);
                    alphaHead = ref Unsafe.Add(ref alphaHead, alpha.Stride);
                }
            }
            else
            {
                Parallel.For(0, batchLength, index =>
                {
                    ref var head = ref Unsafe.Add(ref array.GetHeadRef(),
                        (nint)index * array.Lengths[0]);
                    Span<SingleIndice> indices =
                    stackalloc SingleIndice[array.Rank - 1];
                    for (int i = 0; i < array.Rank - 1; i++)
                    {
                        indices[i] = new SingleIndice(array.Lengths[i + 1],
                            array.Strides[i + 1]);
                    }
                    ApplyWith_Impl<TAction, double>(ref head, alpha[(nint)index],
                        indices);
                });
            }
        }

        public static void ApplyAssign<TAction>(NDArray src, NDArray dest)
            where TAction : struct, IUnaryOperator<double, double>
        {
            if (src.Rank == 0 || dest.Rank == 0)
                return;
            ThrowUtils.ThrowIfArrayNotMatch(src, dest);
            ref var srcHead = ref src.Data[src.Offset];
            ref var destHead = ref dest.Data[dest.Offset];
            Span<DoubleIndice> indices = stackalloc DoubleIndice[src.Rank];
            for (int i = 0; i < src.Rank; i++)
            {
                indices[i] = new DoubleIndice(src.Lengths[i],
                    src.Strides[i], dest.Strides[i]);
            }
            ApplyAssign_Impl<TAction>(ref srcHead, ref destHead, 
                indices);
        }

        private static void ApplyAssign_Impl<TAction>(ref double srcHead, 
            ref double destHead, ReadOnlySpan<DoubleIndice> indices)
            where TAction : struct, IUnaryOperator<double, double>
        {
            if (indices.Length == 0)
                return;
            else if (indices.Length == 1)
            {
                BlasLike.Details.ApplyAssignV_Impl<TAction>( ref srcHead,
                     ref destHead, indices[0]);
            }
            else if (indices.Length == 2)
            {
                BlasLike.Details.ApplyAssignM_Impl<TAction>(ref srcHead, 
                    ref destHead, indices[0], indices[1]);
            }
            else
            {
                var firstIndice = indices[0];
                var subIndices = indices[1..];
                for (var i = 0; i < firstIndice.Length; i++)
                {
                    ApplyAssign_Impl<TAction>(ref srcHead, ref destHead,
                        subIndices);
                    srcHead = ref Unsafe.Add(ref srcHead, firstIndice.AStride);
                    destHead = ref Unsafe.Add(ref destHead, firstIndice.BStride);
                }
            }
        }

        public static void ApplyWithAssign<TAction, TIn>(NDArray src, TIn alpha, NDArray dest)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
        {
            if (src.Rank == 0 || dest.Rank == 0)
                return;
            ThrowUtils.ThrowIfArrayNotMatch(src, dest);
            ref var srcHead = ref src.Data[src.Offset];
            ref var destHead = ref dest.Data[dest.Offset];
            Span<DoubleIndice> indices = stackalloc DoubleIndice[src.Rank];
            for (int i = 0; i < src.Rank; i++)
            {
                indices[i] = new DoubleIndice(src.Lengths[i],
                    src.Strides[i], dest.Strides[i]);
            }
            ApplyWithAssign_Impl<TAction, TIn>(ref srcHead, alpha,
                ref destHead, indices);
        }

        private static void ApplyWithAssign_Impl<TAction, TIn>(ref double srcHead,
            TIn alpha, ref double destHead, ReadOnlySpan<DoubleIndice> indices)
            where TAction : struct, IBinaryOperator<double, TIn, double>
            where TIn : struct
        {
            if (indices.Length == 0)
                return;
            else if (indices.Length == 1)
            {
                BlasLike.Details.ApplyWithAssignV_Impl<TAction, TIn>(ref srcHead,
                    alpha, ref destHead, indices[0]);
            }
            else if (indices.Length == 2)
            {
                BlasLike.Details.ApplyWithAssignM_Impl<TAction, TIn>(ref srcHead,
                    alpha, ref destHead, indices[0], indices[1]);
            }
            else
            {
                var firstIndice = indices[0];
                var subIndices = indices[1..];
                for (var i = 0; i < firstIndice.Length; i++)
                {
                    ApplyWithAssign_Impl<TAction, TIn>(ref srcHead, alpha,
                        ref destHead, subIndices);
                    srcHead = ref Unsafe.Add(ref srcHead, firstIndice.AStride);
                    destHead = ref Unsafe.Add(ref destHead, firstIndice.BStride);
                }
            }
        }

        public static void BatchApplyWithAssign<TAction>(NDArray src, Vector alpha, NDArray dest)
            where TAction : struct, IBinaryOperator<double, double, double>
        {
            if (src.Rank != dest.Rank + 1)
                throw new ArgumentException($"The rank of {nameof(src)} " +
                    $"must be equals to {nameof(dest)}.Rank + 1.");
            var batchLength = src.Lengths[0];
            if (alpha.Length != batchLength)
                throw new ArgumentException($"The length of {nameof(alpha)} " +
                    $"must be equals to {nameof(src)}.Lengths[0].");
            if (!src.Lengths[1..].SequenceEqual(dest.Lengths))
                throw new ArgumentException($"The right lengths of {nameof(src)} " +
                    $"must be equals to {nameof(dest)}.Lengths.");

            ref var srcHead = ref src.GetHeadRef();
            ref var alphaHead = ref alpha.GetHeadRef();
            ref var destHead = ref dest.GetHeadRef();
            Span<DoubleIndice> indices =
                    stackalloc DoubleIndice[src.Rank - 1];
            for (int i = 0; i < src.Rank - 1; i++)
            {
                indices[i] = new DoubleIndice(src.Lengths[i + 1],
                    src.Strides[i + 1], dest.Strides[i]);
            }
            if (src.Rank == 1)
            {
                for (int i = 0; i < batchLength; i++)
                {
                    destHead = TAction.Invoke(alphaHead, srcHead);
                    srcHead = ref Unsafe.Add(ref srcHead, src.Strides[0]);
                    alphaHead = ref Unsafe.Add(ref alphaHead, alpha.Stride);
                }
            }
            else
            {
                for (int i = 0; i < batchLength; i++)
                {
                    ApplyWithAssign_Impl<TAction, double>(ref srcHead,
                        alphaHead, ref destHead, indices);
                    srcHead = ref Unsafe.Add(ref srcHead, src.Strides[0]);
                    alphaHead = ref Unsafe.Add(ref alphaHead, alpha.Stride);
                }
            }
        }

        public static void ApplyTo<TAction>(NDArray src, NDArray dest)
            where TAction : struct, IBinaryOperator<double, double, double>
        {
            if (src.Rank == 0)
                return;
            ThrowUtils.ThrowIfArrayNotMatch(src, dest);
            ref var srcHead = ref src.Data[src.Offset];
            ref var destHead = ref dest.Data[dest.Offset];
            Span<DoubleIndice> indices = stackalloc DoubleIndice[src.Rank];
            for (int i = 0; i < src.Rank; i++)
            {
                indices[i] = new DoubleIndice(src.Lengths[i],
                    src.Strides[i], dest.Strides[i]);
            }
            ApplyTo_Impl<TAction>(ref srcHead, ref destHead, indices);
        }

        private static void ApplyTo_Impl<TAction>(ref double srcHead, 
            ref double destHead, ReadOnlySpan<DoubleIndice> indices)
            where TAction : struct, IBinaryOperator<double, double, double>
        {
            if (indices.Length == 0)
                return;
            else if (indices.Length == 1)
            {
                BlasLike.Details.ApplyToV_Impl<TAction>(ref srcHead, ref destHead,
                    indices[0]);
            }
            else if (indices.Length == 2)
            {
                BlasLike.Details.ApplyToM_Impl<TAction>(ref srcHead, ref destHead,
                    indices[0], indices[1]);
            }
            else
            {
                var firstIndice = indices[0];
                var subIndices = indices[1..];
                for (var i = 0; i < firstIndice.Length; i++)
                {
                    ApplyTo_Impl<TAction>(ref srcHead, ref destHead, subIndices);
                    srcHead = ref Unsafe.Add(ref srcHead, firstIndice.AStride);
                    destHead = ref Unsafe.Add(ref destHead, firstIndice.BStride);
                }
            }
        }

        public static void ApplyToWith<TAction, TIn>(NDArray src, TIn alpha, NDArray dest)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
        {
            if (src.Rank == 0 || dest.Rank == 0)
                return;
            ThrowUtils.ThrowIfArrayNotMatch(src, dest);
            ref var srcHead = ref src.Data[src.Offset];
            ref var destHead = ref dest.Data[dest.Offset];
            Span<DoubleIndice> indices = stackalloc DoubleIndice[src.Rank];
            for (int i = 0; i < src.Rank; i++)
            {
                indices[i] = new DoubleIndice(src.Lengths[i],
                    src.Strides[i], dest.Strides[i]);
            }
            ApplyToWith_Impl<TAction, TIn>(ref srcHead, alpha, 
                ref destHead, indices);
        }

        private static void ApplyToWith_Impl<TAction, TIn>(ref double srcHead, 
            TIn alpha, ref double destHead,ReadOnlySpan<DoubleIndice> indices)
            where TAction : struct, ITernaryOperator<double, TIn, double, double>
            where TIn : struct
        {
            if (indices.Length == 0)
                return;
            else if (indices.Length == 1)
            {
                BlasLike.Details.ApplyToWithV_Impl<TAction, TIn>(ref srcHead,
                    alpha, ref destHead, indices[0]);
            }
            else if (indices.Length == 2)
            {
                BlasLike.Details.ApplyToWithM_Impl<TAction, TIn>(ref srcHead, 
                    alpha, ref destHead, indices[0], indices[1]);
            }
            else
            {
                var firstIndice = indices[0];
                var subIndices = indices[1..];
                for (var i = 0; i < firstIndice.Length; i++)
                {
                    ApplyToWith_Impl<TAction, TIn>(ref srcHead, alpha, 
                        ref destHead, subIndices);
                    srcHead = ref Unsafe.Add(ref srcHead, firstIndice.AStride);
                    destHead = ref Unsafe.Add(ref destHead, firstIndice.BStride);
                }
            }
        }

        public static void BatchApplyToWith<TAction>(NDArray src, Vector alpha, NDArray dest)
            where TAction : struct, ITernaryOperator<double, double, double, double>
        {
            if (src.Rank != dest.Rank + 1)
                throw new ArgumentException($"The rank of {nameof(src)} " +
                    $"must be equals to {nameof(dest)}.Rank + 1.");
            var batchLength = src.Lengths[0];
            if (alpha.Length != batchLength)
                throw new ArgumentException($"The length of {nameof(alpha)} " +
                    $"must be equals to {nameof(src)}.Lengths[0].");
            if (!src.Lengths[1..].SequenceEqual(dest.Lengths))
                throw new ArgumentException($"The right lengths of {nameof(src)} " +
                    $"must be equals to {nameof(dest)}.Lengths.");
            ref var srcHead = ref src.GetHeadRef();
            ref var alphaHead = ref alpha.GetHeadRef();
            ref var destHead = ref dest.GetHeadRef();
            Span<DoubleIndice> indices =
                    stackalloc DoubleIndice[src.Rank - 1];
            for (int i = 0; i < src.Rank - 1; i++)
            {
                indices[i] = new DoubleIndice(src.Lengths[i + 1],
                    src.Strides[i + 1], dest.Strides[i]);
            }
            if (src.Rank == 1)
            {
                for (int i = 0; i < batchLength; i++)
                {
                    destHead = TAction.Invoke(srcHead, alphaHead, destHead);
                    srcHead = ref Unsafe.Add(ref srcHead, src.Strides[0]);
                    alphaHead = ref Unsafe.Add(ref alphaHead, alpha.Stride);
                }
            }
            else
            {
                for (int i = 0; i < batchLength; i++)
                {
                    ApplyToWith_Impl<TAction, double>(ref srcHead,
                        alpha[i], ref destHead, indices);
                    srcHead = ref Unsafe.Add(ref srcHead, src.Strides[0]);
                    alphaHead = ref Unsafe.Add(ref alphaHead, alpha.Stride);
                }
            }
        }
    }
}
