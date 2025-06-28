using SimpleHelpers.Indices;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace SimpleHelpers.LinearAlg
{
    
    [DebuggerTypeProxy(typeof(DVectorSpanDebugView))]
    public class Vector
    {
        public double[] Data { get; }
        public int Offset { get; }
        public nint Length { get; }
        public nint Stride { get; }

        public SingleIndice Indice => new(Length, Stride);


        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector(double[] array)
        {
            ArgumentNullException.ThrowIfNull(array, nameof(array));
            Data = array;
            Offset = 0;
            Length = array.Length;
            Stride = 1;
        }

        
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector(double[] array, int start, nint length, nint step = 1)
        {
            ArgumentNullException.ThrowIfNull(array, nameof(array));
            ArgumentOutOfRangeException.ThrowIfNegative(start, nameof(start));
            ArgumentOutOfRangeException.ThrowIfLessThan(step, 1, nameof(step));
            ArgumentOutOfRangeException.ThrowIfNegative(length, nameof(length));
            ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(array.Length, start + (length - 1) * step);

            Data = array;
            Offset = start;
            Length = length;
            Stride = step;
        }

        public static Vector Create(nint length, bool uninited = false)
        {
            ArgumentOutOfRangeException.ThrowIfGreaterThan(length, Array.MaxLength, nameof(length));
            double[] data = uninited ? GC.AllocateUninitializedArray<double>((int)length)
                : new double[length];
            return new(data);
        }

        public static Vector Create(nint length, double val)
        {
            ArgumentOutOfRangeException.ThrowIfGreaterThan(length, Array.MaxLength, nameof(length));
            double[] data = GC.AllocateUninitializedArray<double>((int)length);
            data.AsSpan().Fill(val);
            return new(data);
        }

        public ref double this[nint index] => ref At(index);

        public Vector this[NRange range]
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get
            {
                (var start, var length) = range.GetOffsetAndLength(Length);
                return Slice(start, length);
            }
        }

        public bool IsEmpty => Length == 0;

        public override bool Equals(object? obj) =>
                throw new NotSupportedException();

        public override int GetHashCode() =>
            throw new NotSupportedException();

        
        public static Span<double> Empty => default;

        
        public VectorEnumerator GetEnumerator() => new(this);

        public Vector MakeContinous(bool forceCopy = false) 
        {
            if(!forceCopy && Stride == 1)
                return this;
            return Clone();
        }

        internal Vector Clone()
        {
            if (IsEmpty)
                return new Vector([]);
            var target = new Vector(new double[Length], 0, Length, Stride);
            BlasLike.Copy(this, target);
            return target;
        }

        internal ref double GetHeadRef()
        {
            if (IsEmpty)
                throw new InvalidOperationException("Matrix is empty.");
            return ref Data[Offset];
        }

        internal Span<double> GetSpan()
        {
            if (IsEmpty)
                throw new InvalidOperationException("Matrix is empty.");
            return Data.AsSpan(Offset);
        }

        public ref struct VectorEnumerator
        {
            private readonly Vector _span;
            
            private nint _index;

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            internal VectorEnumerator(Vector span)
            {
                _span = span;
                _index = -1;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public bool MoveNext()
            {
                var index = _index + 1;
                if (index < _span.Length)
                {
                    _index = index;
                    return true;
                }

                return false;
            }

            public readonly ref double Current
            {
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                get => ref _span[_index];
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe void Clear()
        {
            if (Length == 0)
            {
                return;
            }
            ref double current = ref Data[Offset];
            if (Stride == 1)
            {
                for (nint i = 0; i < Length; i += int.MaxValue / 2)
                {
                    int length = (int)Math.Min(int.MaxValue / 2, Length - i);
                    var span = MemoryMarshal.CreateSpan(ref current, length);
                    span.Clear();
                    current = ref Unsafe.Add(ref current, length);
                }
            }
            else
            {
                for (nint i = 0; i < Length; i++)
                {
                    current = default;
                    current = ref Unsafe.Add(ref current, Stride);
                }
            }
            return;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe void Fill(double value)
        {
            if (Length == 0)
                return;
            if (value == default)
            {
                Clear();
                return;
            }

            ref double current = ref Data[Offset];
            if (Stride == 1)
            {
                for (nint i = 0; i < Length; i += int.MaxValue / 2)
                {
                    int length = (int)Math.Min(int.MaxValue / 2, Length - i);
                    var span = MemoryMarshal.CreateSpan(ref current, length);
                    span.Fill(value);
                    current = ref Unsafe.Add(ref current, length);
                }
            }
            else
            {
                for (nint i = 0; i < Length; i++)
                {
                    current = value;
                    current = ref Unsafe.Add(ref current, Stride);
                }
            }
            return;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector Slice(nint start)
        {
            if (start > Length)
                throw new ArgumentOutOfRangeException(nameof(start), start, "Error: start index should be less than length.");

            return new Vector(Data, (int)(Offset + start * Stride), Length - start, Stride);
        }

        public void CopyTo(Vector other)
        {
            ArgumentNullException.ThrowIfNull(other, nameof(other));
            ArgumentOutOfRangeException.
                ThrowIfLessThan(other.Length, Length, nameof(other));
            ref double current = ref Data[Offset];
            ref double otherCurrent = ref other.Data[other.Offset];
            for (nint i = 0; i < Length; i++)
            {
                otherCurrent = current;
                current = ref Unsafe.Add(ref current, Stride);
                otherCurrent = ref Unsafe.Add(ref otherCurrent, other.Stride);
            }
        }

        public ref double At(nint index)
        {
            if (index < 0 || index >= Length)
                throw new IndexOutOfRangeException(nameof(index));
            return ref AtUncheck(index);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal ref double AtUncheck(nint index)
        {
            return ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(Data), Offset + index * Stride);
        }


        
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector Slice(nint start, nint length)
        {
            if (start > Length)
                throw new ArgumentOutOfRangeException(nameof(start), start, "Error: start index should be less than length.");
            if (length < 0 || start + length > Length)
                throw new ArgumentOutOfRangeException(nameof(length), length, "Error: length should be greater than 0 and end of span should be in range.");
            return SliceUncheck(start, length);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal Vector SliceUncheck(nint start, nint length)
        {
            return new Vector(Data, (int)(Offset + start * Stride), length, Stride);
        }

        public void FlattenTo(Span<double> target)
        {
            if (target.Length < Length)
                throw new ArgumentException("Error: target length should be greater than source length.", nameof(target));
            ref double current = ref Data[Offset];
            int length = (int)Length;
            if (Stride == 1)
            {
                MemoryMarshal.CreateSpan(ref current, length).CopyTo(target);
                return;
            }
            ref double targetCurrent = ref MemoryMarshal.GetReference(target);
            for (int i = 0; i < length; i++)
            {
                targetCurrent = current;
                current = ref Unsafe.Add(ref current, Stride);
                targetCurrent = ref Unsafe.Add(ref targetCurrent, 1);
            }
        }

        public static double operator *(Vector left, Vector right)
            => BlasLike.Dot(left, right);
    }

    public readonly ref struct DVectorSpanDebugView
    {
        public readonly object Items;
        public DVectorSpanDebugView(Vector span)
        {
            Length = span.Length;
            Stride = span.Stride;

            if (span.IsEmpty || TooLong)
            {
                Items = Array.Empty<double>();
            }
            else
            {
                var data = new double[span.Length];
                span.FlattenTo(data);
                Items = data;
            }
        }

        public long Length { get; }
        public long Stride { get; }

        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public bool TooLong => Length > TooLongLimit;

        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        public const int TooLongLimit = int.MaxValue / 2;
    }
}
