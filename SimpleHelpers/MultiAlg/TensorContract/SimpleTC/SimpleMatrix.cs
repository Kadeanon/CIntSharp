using CommunityToolkit.HighPerformance;
using SimpleHelpers.Indices;
using System.Runtime.CompilerServices;

namespace SimpleHelpers.MultiAlg.TensorContract.SimpleTC
{
    internal readonly struct SimpleMatrix // ref struct means we gave up
    {                                         // the ability to multi-threading.
        readonly double[] data;             // dont use it for large system.
        readonly nint offset;
        public readonly nint[] rowScatters;
        public readonly nint[] colScatters;

        public SimpleMatrix(NDArray array,
            SingleIndice[] rowIndices, SingleIndice[] colIndices)
        {
            Span<nint> zeros = stackalloc nint[array.Rank];
            this.data = array.Data;
            offset = array.Offset;
            rowScatters = BuildScatters(rowIndices);
            colScatters = BuildScatters(colIndices);
        }

        public nint LocalI(nint row)
            => rowScatters[row];

        public nint LocalJ(nint col)
            => colScatters[col];

        public ref double Local(nint row, nint col)
        {
            nint rowOffset = LocalI(row);
            nint colOffset = LocalJ(col);
            return ref data.DangerousGetReferenceAt
                ((int)(offset + rowOffset + colOffset));
        }

        public ref double this[nint row, nint col]
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => ref Local(row, col);
        }

        public nint RowLength => rowScatters.Length;

        public nint ColLength => colScatters.Length;

        public static nint[] BuildScatters(SingleIndice[] indices)
        {

            nint size = 1;
            for (int i = 0; i < indices.Length; i++)
            {
                size *= indices[i].Length;
            }
            var scatters = new nint[size];

            for (nint i = 0; i < size; i++)
            {
                var index = i;
                nint offset = 0;
                for (nint j = indices.Length - 1; j >= 0; j--)
                {
                    var indice = indices[j];
                    (index, nint currentIndex) = Math.DivRem(index, indice.Length);
                    offset += currentIndex * indice.Stride;
                }
                scatters[i] = offset;
            }

            return scatters;
        }
    }
}
