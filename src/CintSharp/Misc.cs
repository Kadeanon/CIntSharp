using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp
{
    internal static class Misc
    {
        public static void MALLOC_INSTACK<TTo>(scoped ref Span<double> buffer, out Span<TTo> split, int n)
            where TTo : unmanaged
        {
            Span<byte> bytesBuffer = MemoryMarshal.AsBytes(buffer);
            int spaceToUse = n * Unsafe.SizeOf<TTo>();
            int total = bytesBuffer.Length;
            if(spaceToUse > total)
            {
                throw new ArgumentException("Not enough space in buffer");
            }
            split = MemoryMarshal.Cast<byte, TTo>(bytesBuffer[..spaceToUse]);
            bytesBuffer = bytesBuffer[spaceToUse..];
            buffer = MemoryMarshal.Cast<byte, double>(bytesBuffer);
        }
    }
}
