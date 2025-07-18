using SimpleHelpers.Utilities.Pools;
using System.ComponentModel;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using Zyl.ExSpans;

namespace SimpleHelpers.LinearAlg
{
    /// <summary>
    /// <see cref="VectorSpan"/> is a type that represents a contiguous region of arbitrary memory with step, similar to <see cref="Span{double}"/>.
    /// </summary>
    public readonly ref struct VectorSpan : IEquatable<VectorSpan>
    {
        private readonly ref double _reference;
        private readonly nint _length;
        private readonly nint _step;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public VectorSpan(double[]? array)
        {
            if (array == null)
            {
                this = default;
                return; // returns default
            }
            if (!typeof(double).IsValueType && array.GetType() != typeof(double[]))
                throw new ArgumentException("Error: array is not a double[] object.", nameof(array));

            _reference = ref MemoryMarshal.GetArrayDataReference(array);
            _length = array.Length;
            _step = 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public VectorSpan(double[] array, nint start, nint length, nint step = 1)
        {
            if (array == null)
            {
                throw new ArgumentException("Error: array should be not null.", nameof(array));
            }
            if (!typeof(double).IsValueType && array.GetType() != typeof(double[]))
                throw new ArgumentException("Error: array is not a double[] object.", nameof(array));

            if (start > array.Length)
                throw new ArgumentOutOfRangeException(nameof(start), start, "Error: end of span is out of range.");
            var totalSize = length * step;
            if (start + totalSize > array.Length)
                throw new ArgumentOutOfRangeException(nameof(length), length, "Error: end of span is out of range.");

            _reference = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(array), (nint)(uint)start);
            _length = length;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe VectorSpan(double* pointer, nint length, nint step)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(length, nameof(length));

            _reference = ref Unsafe.AsRef<double>(pointer);
            _length = length;
            _step = step;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public VectorSpan(ref double reference)
        {
            _reference = ref reference;
            _length = 1;
            _step = 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal VectorSpan(ref double reference, nint length, nint step = 1)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(length, nameof(length));

            _reference = ref reference;
            _length = length;
            _step = step;
        }

        public unsafe ref double this[nint index]
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get
            {
                if ((ulong)index >= (ulong)_length)
                    throw new IndexOutOfRangeException(nameof(index));
                return ref Unsafe.Add(ref _reference, index * _step);
            }
        }

        public nint Length => _length;

        public nint Stride => _step;

        public bool IsEmpty
        {
            get => _length == 0;
        }

        public static bool operator !=(VectorSpan left, VectorSpan right) => !(left == right);

        public override bool Equals(object? obj) =>
                throw new NotSupportedException();

        public override int GetHashCode() =>
            throw new NotSupportedException();

        public static implicit operator VectorSpan(double[] array) => new VectorSpan(array);

        public static Span<double> Empty => default;

        public NStepSpanEnumerator GetEnumerator() => new NStepSpanEnumerator(this);

        public ref struct NStepSpanEnumerator
        {
            /// <summary>The span being enumerated.</summary>
            private readonly VectorSpan _span;
            /// <summary>The next index to yield.</summary>
            private nint _index;

            /// <summary>Initialize the enumerator.</summary>
            /// <param name="span">The span to enumerate.</param>
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            internal NStepSpanEnumerator(VectorSpan span)
            {
                _span = span;
                _index = -1;
            }

            /// <summary>Advances the enumerator to the next element of the span.</summary>
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

            /// <summary>Gets the element at the current position of the enumerator.</summary>
            public ref double Current
            {
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                get => ref _span[_index];
            }
        }

        [Browsable(false)]
        public ref double GetPinnableReference()
        {
            // Ensure that the native code has just one forward branch that is predicted-not-taken.
            ref double ret = ref Unsafe.NullRef<double>();
            if (_length != 0) ret = ref _reference;
            return ref ret;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe void Clear()
        {
            if (_step == 1)
            {
                if (_length == 0)
                {
                    return;
                }
                var span = ExMemoryMarshal.CreateExSpan
                    (ref _reference, _length);
                span.Clear();
                return;
            }

            nint index = 0;
            ref double current = ref _reference;
            while (index < _length)
            {
                current = default;
                current = ref Unsafe.Add(ref current, _step);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe void Fill(double value)
        {
            if (_step == 1)
            {
                if (_length == 0)
                {
                    return;
                }
                var span = ExMemoryMarshal.CreateExSpan
                    (ref _reference, _length);
                span.Fill(value);
                return;
            }

            nint index = 0;
            ref double current = ref _reference;
            while (index < _length)
            {
                current = value;
                current = ref Unsafe.Add(ref current, _step);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe void CopyTo(ExSpan<double> destination)
        {
            if (_length == 0)
            {
                return;
            }
            if (_step == 1)
            {
                var span = ExMemoryMarshal.CreateExSpan
                    (ref _reference, _length);
                span.CopyTo(destination);
                return;
            }

            for (nint i = 0; i < _length; i++)
            {
                destination[i] = this[i];
            }
        }

        public unsafe void CopyTo(VectorSpan destination)
        {
            if (_length == 0)
            {
                return;
            }
            if (destination._step == 1)
            {
                var span = ExMemoryMarshal.CreateExSpan
                    (ref destination._reference, destination._length);
                CopyTo(span);
                return;
            }

            ref double src = ref _reference;
            ref double dst = ref destination._reference;
            nint length = Math.Min(_length, destination._length);
            for (nint i = 0; i < length; i++)
            {
                destination[i] = this[i];
                dst = ref Unsafe.Add(ref dst, destination._step);
                src = ref Unsafe.Add(ref src, _step);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe bool CopyToSafe(VectorSpan destination)
        {
            if (destination.Length < _length)
            {
                var span = Slice(0, destination.Length);
                span.CopyTo(destination);
                return false;
            }

            CopyTo(destination);
            return true;
        }

        public static bool operator ==(in VectorSpan left, in VectorSpan right)
        {
            if (left._length != right._length)
                return false;

            ref double leftRef = ref left._reference;
            ref double rightRef = ref right._reference;
            nint unrollSize = 4;
            nint length = left._length;
            nint i = 0;
            if (length >= 2 * unrollSize)
            {
                nint leftStride = left._step;
                nint leftStride2 = leftStride * 2;
                nint leftStride3 = leftStride + leftStride2;
                nint leftStride4 = leftStride2 + leftStride2;
                nint rightStride = right._step;
                nint rightStride2 = rightStride * 2;
                nint rightStride3 = rightStride + rightStride2;
                nint rightStride4 = rightStride2 + rightStride2;
                for (; i < length - unrollSize; i += unrollSize)
                {
                    bool neq = false;
                    neq |= leftRef == rightRef;
                    neq |= Unsafe.Add(ref leftRef, leftStride)
                        == Unsafe.Add(ref rightRef, rightStride);
                    neq |= Unsafe.Add(ref leftRef, leftStride2)
                        == Unsafe.Add(ref rightRef, rightStride2);
                    neq |= Unsafe.Add(ref leftRef, leftStride3)
                        == Unsafe.Add(ref rightRef, rightStride3);
                    if (neq)
                        return false;
                    leftRef = ref Unsafe.Add(ref leftRef, leftStride4);
                    rightRef = ref Unsafe.Add(ref rightRef, rightStride4);
                }
            }
            for (; i < length; i++)
            {
                if (leftRef != rightRef)
                    return false;
                leftRef = ref Unsafe.Add(ref leftRef, left._step);
                rightRef = ref Unsafe.Add(ref rightRef, right._step);
            }
            return true;
        }

        public override string ToString()
        {
            var sb = StringBuilderPool.Shared.Get();
            ToString(sb, 6, 4);
            var result = sb.ToString();
            StringBuilderPool.Shared.Return(sb);
            return $"VectorSpan[{_length}]\n[{result}]";
        }

        /// <summary>
        /// Returns a string representation of the elements in the span.
        /// </summary>
        /// <param name="sb"> The <see cref="StringBuilder"/> object to write the string</param>
        /// <param name="start">If the span is too long, determines how many element to write before '...'.</param>
        /// <param name="end">If the span is too long, determines how many element to write after '...'.</param>
        public void ToString(StringBuilder sb, nint start, nint end, string? format = null)
        {
            if (_length == 0)
            {
                return;
            }
            ref double current = ref _reference;

            if (format is not null)
            {
                if (_length <= start + end)
                {
                    sb.Append(current.ToString(format));
                    nint i = 1;
                    for (; i < _length; i++)
                    {
                        current = ref Unsafe.Add(ref current, _step);
                        sb.Append(", ");
                        sb.Append(current.ToString(format));
                    }
                }
                else
                {
                    for (nint i = 0; i < start; i++)
                    {
                        if (i > 0) sb.Append(", ");
                        sb.Append(this[i].ToString(format));
                    }
                    sb.Append(", ... ");
                    for (nint i = _length - end; i < _length; i++)
                    {
                        if (i > 0) sb.Append(", ");
                        sb.Append(this[i].ToString(format));
                    }
                }
            }
            else
            {
                if (_length <= start + end)
                {
                    sb.Append(current);
                    nint i = 1;
                    for (; i < _length; i++)
                    {
                        current = ref Unsafe.Add(ref current, _step);
                        sb.Append(", ");
                        sb.Append(current);
                    }
                }
                else
                {
                    for (nint i = 0; i < start; i++)
                    {
                        if (i > 0) sb.Append(", ");
                        sb.Append(this[i]);
                    }
                    sb.Append(", ... ");
                    for (nint i = _length - end; i < _length; i++)
                    {
                        if (i > 0) sb.Append(", ");
                        sb.Append(this[i]);
                    }
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public VectorSpan Slice(nint start)
        {
            if (start > _length)
                throw new ArgumentOutOfRangeException(nameof(start), start, "Error: start index should be less than length.");

            return new VectorSpan(ref this[start], _length - start, _step);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public VectorSpan Slice(nint start, nint length)
        {
            if (start > _length)
                throw new ArgumentOutOfRangeException(nameof(start), start, "Error: start index should be less than length.");
            if (length < 0 || start + length > _length)
                throw new ArgumentOutOfRangeException(nameof(length), length, "Error: length should be greater than 0 and end of span should be in range.");
            return new VectorSpan(ref this[start], length, _step);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double[] ToArray()
        {
            if (_length == 0)
                return Array.Empty<double>();

            var destination = new double[_length];
            CopyTo(new ExSpan<double>(destination));
            return destination;
        }

        internal unsafe ExSpan<double> AsExSpan()
        {
            if (_length == 0)
                return ExSpan<double>.Empty;
            if (_step != 1)
                throw new NotSupportedException("Error: step should be 1.");

            return ExMemoryMarshal.CreateExSpan(ref _reference, _length);
        }

        internal ExSpan<double> AsExSpan(nint start, int length)
        {
            if (start > _length)
                throw new ArgumentOutOfRangeException(nameof(start), start, "Error: start index should be less than length.");
            if (length < 0 || start + length > _length)
                throw new ArgumentOutOfRangeException(nameof(length), length, "Error: length should be greater than 0 and end of span should be in range.");
            if (_step != 1)
                throw new NotSupportedException("Error: step should be 1.");

            if (length == 0)
                return ExSpan<double>.Empty;
            if (_step != 1)
                throw new NotSupportedException("Error: step should be 1.");

            ref double target = ref Unsafe.Add(ref _reference, start);
            return ExMemoryMarshal.CreateExSpan(ref target, length);
        }

        internal ReadOnlyExSpan<double> AsReadOnlyExSpan()
        {
            if (_length == 0)
                return ReadOnlyExSpan<double>.Empty;
            if (_step != 1)
                throw new NotSupportedException("Error: step should be 1.");

            return ExMemoryMarshal.CreateReadOnlyExSpan
                (ref _reference, _length);
        }

        internal ReadOnlyExSpan<double> AsReadOnlyExSpan(nint start, int length)
        {
            if (start > _length)
                throw new ArgumentOutOfRangeException(nameof(start), start, "Error: start index should be less than length.");
            if (length < 0 || start + length > _length)
                throw new ArgumentOutOfRangeException(nameof(length), length, "Error: length should be greater than 0 and end of span should be in range.");
            if (_step != 1)
                throw new NotSupportedException("Error: step should be 1.");

            if (length == 0)
                return ExSpan<double>.Empty;
            if (_step != 1)
                throw new NotSupportedException("Error: step should be 1.");

            ref double target = ref Unsafe.Add(ref _reference, start);
            return ExMemoryMarshal.CreateReadOnlyExSpan(ref target, length);
        }


        internal unsafe double* Address
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get
            {
                return (double*)Unsafe.AsPointer(ref _reference);
            }
        }

        public bool Equals(VectorSpan other)
            => _length == other._length &&
               _step == other._step &&
               Unsafe.AreSame(ref _reference, ref other._reference);

        public static implicit operator VectorSpan(Span<double> span) => new(ref span[0], span.Length);

        public static implicit operator VectorSpan(ExSpan<double> span) => new(ref span[0], span.Length);
    }
}
