using SimpleHelpers.MultiAlg;
using SimpleHelpers.Utilities.Pools;
using System.Buffers;
using System.Diagnostics;
using System.Numerics;

namespace SimpleHelpers.Utilities
{
    internal static class ThrowUtils
    {

        [DebuggerHidden]
        internal static void ThrowIfNotInRange_LeftClosedRightClosed<T>(T value, T min, T max, string name)
            where T : IComparable<T>
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(value, min, name);
            ArgumentOutOfRangeException.ThrowIfGreaterThan(value, max, name);
        }

        [DebuggerHidden]
        internal static void ThrowIfNotInRange_LeftClosedRightOpen<T>(T value, T min, T max, string name)
            where T : IComparable<T>
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(value, min, name);
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(value, max, name);
        }

        [DebuggerHidden]
        internal static void ThrowIfNotInRange_LeftOpenRightClosed<T>(T value, T min, T max, string name)
            where T : IComparable<T>
        {
            ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(value, min, name);
            ArgumentOutOfRangeException.ThrowIfGreaterThan(value, max, name);
        }

        [DebuggerHidden]
        internal static void ThrowIfNotInRange_LeftOpenRightOpen<T>(T value, T min, T max, string name)
            where T : IComparable<T>
        {
            ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(value, min, name);
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(value, max, name);
        }

        [DebuggerHidden]
        internal static void ThrowIfNotEqual<T>(T value, T expected, string name)
            where T : IEquatable<T>
        {
            if (!value.Equals(expected))
            {
                throw new ArgumentException($"The {name} is not equal to the expected value {expected}.");
            }
        }

        [DebuggerHidden]
        internal static void ThrowIfRankNotMatch<T>(ReadOnlySpan<T> value, int rank, string name)
        {
            if (value.Length != rank)
            {
                throw new ArgumentException($"The rank of the {name} does not match the rank of the NDArray.");
            }
        }

        /// <summary>
        /// Throws an exception if the value is not in the range of [0, max).
        /// If <paramref name="value"/> is less than 0, it will be added to <paramref name="max"/> to check if it is in the range.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="value"></param>
        /// <param name="max"></param>
        /// <param name="name"></param>
        /// <returns></returns>
        [DebuggerHidden]
        internal static void ThrowIfNotInRange_CheckNegative<T>(ref T value, T max, string name)
            where T : INumber<T>
        {
            if (T.IsNegative(value))
            {
                value += max;
            }
            ArgumentOutOfRangeException.ThrowIfNegative(value, name);
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(value, max, name);
        }

        [DebuggerHidden]
        internal static nint ThrowIfNotInRange_CheckNIndex(NIndex range, nint max, string name)
        {
            nint offset = range.GetOffset(max);
            ArgumentOutOfRangeException.ThrowIfNegative(offset, name);
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(offset, max, name);
            return offset;
        }

        [DebuggerHidden]
        internal static (nint offset, nint length) ThrowIfNotInRange_CheckNRange(NRange range, nint max, string name)
        {
            var (offset, length) = range.GetOffsetAndLength(max);
            ArgumentOutOfRangeException.ThrowIfNegative(offset, name);
            ArgumentOutOfRangeException.ThrowIfNegative(length, name);
            ArgumentOutOfRangeException.ThrowIfGreaterThan(offset + length, max, name);
            return (offset, length);
        }

        [DebuggerHidden]
        internal static void ThrowIfArrayRankNotEqual<T>(int rank, NDArray array)
            where T : unmanaged, INumberBase<T>
        {
            if (array.Rank != rank)
            {
                throw new ArgumentException($"The rank of the array does not match the rank of the NDArray.");
            }
        }

        [DebuggerHidden]
        internal static void ThrowIfArrayLengthNotEqual<T>(nint length, NDArray array)
            where T : unmanaged, INumberBase<T>
        {
            if (array.Size != length)
            {
                throw new ArgumentException($"The total length of the array does not match the rank of the NDArray.");
            }
        }

        [DebuggerHidden]
        internal static void ThrowIfArrayNotMatch(NDArray left, NDArray right)
        {
            bool areEqual = true;
            if (left.Rank != right.Rank)
            {
                areEqual= false;
            }
            if (areEqual) 
            { 
                for (int i = 0; i < left.Lengths.Length; i++)
                {
                    if (left.Lengths[i] != right.Lengths[i])
                    {
                        areEqual = false;
                        break;
                    }
                }
            }
            if (!areEqual)
            {
                int leftRank = left.Rank;
                int rightRank = right.Rank;
                using var handle = StringBuilderPool.Borrow(out var sb);
                sb.Append("The NDArray lengths do not match: left is (");
                for (int i = 0; i < leftRank; i++)
                {
                    sb.Append(left.Lengths[i]);
                    if (i < leftRank - 1)
                    {
                        sb.Append(", ");
                    }
                }
                sb.Append(") and right is (");
                for (int i = 0; i < rightRank; i++)
                {
                    sb.Append(right.Lengths[i]);
                    if (i < rightRank - 1)
                    {
                        sb.Append(", ");
                    }
                }
                sb.Append(").");
                throw new ArgumentException(sb.ToString());
            }
        }
    }
}
