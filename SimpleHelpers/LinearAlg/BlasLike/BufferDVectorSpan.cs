using NetFabric.Numerics.Tensors.Operators;
using SimpleHelpers.Indices;
using System.Buffers;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace SimpleHelpers.LinearAlg
{
    public unsafe ref struct BufferDVectorSpan : IDisposable
    {
        public ref double bufferHead;
        public ref double oldHead;
        public double[]? managedBuffer;
        public nint oldLength;
        public nint oldStride;
        private bool bufferAllocated;
        private readonly bool shouldCopyBack;

        private readonly bool Managed => oldLength <= int.MaxValue / 2;

        public BufferDVectorSpan(ref double head, nint length, nint stride, bool shouldCopyBack = false, bool forceCopy = false)
        {
            bufferHead = ref head;
            oldHead = ref head;
            oldStride = stride;
            oldLength = length;
            if (length > 0 && stride > 1 || forceCopy)
            {
                this.shouldCopyBack = shouldCopyBack;
                if (Managed)
                {
                    int bufferLength = (int)length;
                    managedBuffer = ArrayPool<double>.Shared.Rent(bufferLength);
                    managedBuffer.AsSpan(0, bufferLength).Clear();
                    bufferHead = ref managedBuffer[0];
                }
                else
                {
                    bufferHead = ref Unsafe.AsRef<double>(
                    NativeMemory.AlignedAlloc(
                    (nuint)(length * sizeof(double)),
                    (nuint)Unsafe.SizeOf<Vector256<double>>()));
                }
                bufferAllocated = true;
                DoubleIndice indice = new(length, stride, 1);
                BlasLike.Details.ApplyAssignV_Impl<IdentityOperator<double>>
                (ref head, ref bufferHead, indice);
            }
        }

        public BufferDVectorSpan(ref double head, nint length, nint stride, double scale, bool shouldCopyBack = false)
        {
            bufferHead = ref head;
            oldHead = ref head;
            oldStride = stride;
            oldLength = length;
            this.shouldCopyBack = shouldCopyBack;
            if (length == 0)
            {
                shouldCopyBack = false;
                bufferAllocated = false;
                return;
            }

            if (Managed)
            {
                int bufferLength = (int)length;
                managedBuffer = ArrayPool<double>.Shared.Rent((int)length);
                managedBuffer.AsSpan(0, bufferLength).Clear();
                bufferHead = ref managedBuffer[0];
            }
            else
            {
                bufferHead = ref Unsafe.AsRef<double>(
                NativeMemory.AlignedAlloc(
                (nuint)(length * sizeof(double)),
                (nuint)Unsafe.SizeOf<Vector256<double>>()));
            }
            bufferAllocated = true;
            DoubleIndice indice = new(length, 1, stride);
            if (scale != 0.0)
                BlasLike.Details.ApplyWithAssignV_Impl
                    <MultiplyOperator<double>, double>(
                    ref head, scale, ref bufferHead, indice);
        }

        public void Dispose()
        {
            if (bufferAllocated)
            {
                if (shouldCopyBack)
                {
                    DoubleIndice indice = new(oldLength, 1, oldStride);
                    BlasLike.Details.ApplyAssignV_Impl<IdentityOperator<double>>
                        (ref bufferHead, ref oldHead, indice);
                }
                if (Managed)
                {
                    int bufferLength = (int)oldLength;
                    managedBuffer.AsSpan(0, bufferLength).Clear();
                    ArrayPool<double>.Shared.Return(managedBuffer!);
                }
                else
                    NativeMemory.AlignedFree(Unsafe.AsPointer(ref bufferHead));
                bufferAllocated = false;
            }
        }
    }
}
