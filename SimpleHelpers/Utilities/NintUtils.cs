using NetFabric.Numerics.Tensors.Operators;
using NetFabric.Numerics.Tensors;

namespace SimpleHelpers.Utilities
{
    public static class NintUtils
    {
        public static bool Is64Bit => Environment.Is64BitProcess;

        //public static nint Product(this IEnumerable<nint> values) => values.Aggregate((nint)1, (a, b) => a * b);

        public static nint Product(this ReadOnlySpan<nint> values)
        => TensorOperations.Product(values) ?? 0;

        public static nint Product(this nint[] values)
        => TensorOperations.Product<nint>(values) ?? 0;

        public static nint Product(this Span<nint> values)
        => TensorOperations.Product<nint>(values) ?? 0;

        public static nint Dot(ReadOnlySpan<nint> left, ReadOnlySpan<nint> right)
        => Tensor.AggregateNumber<nint, MultiplyOperator<nint>, SumOperator<nint>>(left, right);

        public static nint CalculateTotalLength(ReadOnlySpan<nint> lengths, ReadOnlySpan<nint> strides)
        {
            if (lengths.Length != strides.Length)
            {
                throw new ArgumentException("The length of the two spans must be equal.");
            }
            nint sum = 0;
            int index = 0;
            foreach (var value in lengths)
            {
                sum += (value - 1) * strides[index];
                index++;
            }
            return sum;
        }


        public static int IncrementIndexLeft(int curIndex, Span<nint> curIndexes, scoped ReadOnlySpan<nint> length)
        {
            if (curIndex < 0)
                return int.MinValue;
            curIndexes[curIndex] += 1;

            nint max = length[curIndex];
            if (curIndexes[curIndex] < max)
                return 1;
            curIndexes[curIndex] = 0;
            return IncrementIndexLeft(curIndex - 1, curIndexes, length) + 1;
        }

        public static int AddIndexLeft(int curIndex, Span<nint> curIndexes, scoped ReadOnlySpan<nint> length, nint add)
        {
            if (curIndex < 0)
                return 0;
            if (add == 0)
            {
                return 0;
            }
            if (add == 1)
            {
                return IncrementIndexLeft(curIndex, curIndexes, length);
            }
            ref nint current = ref curIndexes[curIndex];
            current += add;
            nint max = length[curIndex];
            if (current < max)
                return 1;
            (nint quotient, nint remainder) = Math.DivRem(curIndexes[curIndex], max);
            current = remainder;
            return AddIndexLeft(curIndex - 1, curIndexes, length, quotient) + 1;
        }
    }
}
