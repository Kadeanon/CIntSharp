using MKLNET;
using NetFabric.Numerics.Tensors;
using NetFabric.Numerics.Tensors.Operators;
using SimpleHelpers.Indices;
using SimpleHelpers.LinearAlg.Dec;
using SimpleHelpers.LinearAlg.Decs;
using System;
using System.Buffers;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Text;
using System.Threading.Tasks;
using Zyl.ExSpans;

namespace SimpleHelpers.LinearAlg
{
    /// <summary>
    /// A span of a matrix with a reference to the first element.
    /// </summary>
    /// <typeparam name="double">The type of the elements in the matrix.</typeparam>
    [DebuggerTypeProxy(typeof(DMatrixSpanDebugView))]
    public class Matrix : IEnumerable<double>
    {
        internal double[] Data { get; }
        internal int Offset { get; set; }
        public nint Rows { get; internal set; }
        public nint Cols { get; internal set; }
        public nint RowStride { get; internal set; }
        public nint ColStride { get; internal set; }

        public SingleIndice RowIndice => new SingleIndice(Rows, RowStride);

        public SingleIndice ColIndice => new SingleIndice(Cols, ColStride);

        public bool IsEmpty => Rows == 0 || Cols == 0;

        public nint TotalSize => Cols * Rows;

        public nint Length => Rows * Cols;

        #region Constructors
        public Matrix(double[] array)
        {
            ArgumentNullException.ThrowIfNull(array, nameof(array));
            Data = array;
            Offset = 0;
            nint length = array.Length;
            Rows = length;
            Cols = 1;
            RowStride = 1;
            ColStride = 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix(double[] array, int offset)
        {
            ArgumentNullException.ThrowIfNull(array, nameof(array));
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual
                (offset, array.Length, nameof(offset));
            Data = array;
            Offset = offset;
            nint length = array.Length - offset;
            Rows = length;
            Cols = 1;
            RowStride = 1;
            ColStride = 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix(double[] array, int rows, nint cols)
        {
            ArgumentNullException.ThrowIfNull(array, nameof(array));
            ArgumentOutOfRangeException.ThrowIfNegative(rows, nameof(rows));
            ArgumentOutOfRangeException.ThrowIfNegative(cols, nameof(cols));
            nint length = rows * cols;
            ArgumentOutOfRangeException.ThrowIfLessThan(array.Length, length, nameof(array));
            Data = array;
            Offset = 0;
            Rows = rows;
            Cols = cols;
            RowStride = cols;
            ColStride = 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix(double[] array, int offset, nint rows, nint cols)
        {
            ArgumentNullException.ThrowIfNull(array, nameof(array));
            ArgumentOutOfRangeException.ThrowIfNegative(offset, nameof(offset));
            ArgumentOutOfRangeException.ThrowIfNegative(rows, nameof(rows));
            ArgumentOutOfRangeException.ThrowIfNegative(cols, nameof(cols));
            nint length = rows * cols;
            ArgumentOutOfRangeException.ThrowIfLessThan(array.Length, offset + length, 
                nameof(array));
            Data = array;
            Offset = offset;
            Rows = rows;
            Cols = cols;
            RowStride = cols;
            ColStride = 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix(double[] array, int offset, nint rows, nint cols, nint rowStride)
        {
            ArgumentNullException.ThrowIfNull(array, nameof(array));
            ArgumentOutOfRangeException.ThrowIfNegative(offset, nameof(offset));
            ArgumentOutOfRangeException.ThrowIfNegative(rows, nameof(rows));
            ArgumentOutOfRangeException.ThrowIfNegative(cols, nameof(cols));
            nint length = (rows - 1) * rowStride + cols;
            ArgumentOutOfRangeException.ThrowIfLessThan(array.Length, offset + length,
                nameof(array));
            Data = array;
            Offset = offset;
            Rows = rows;
            Cols = cols;
            RowStride = rowStride;
            ColStride = 1;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix(double[] array, int offset, nint rows, nint cols,
            nint rowStride, nint colStride)
        {
            ArgumentNullException.ThrowIfNull(array, nameof(array));
            ArgumentOutOfRangeException.ThrowIfNegative(offset, nameof(offset));
            ArgumentOutOfRangeException.ThrowIfNegative(rows, nameof(rows));
            ArgumentOutOfRangeException.ThrowIfNegative(cols, nameof(cols));
            nint last = offset + (rows - 1) * rowStride + (cols - 1) * colStride + 1;
            ArgumentOutOfRangeException.ThrowIfLessThan(array.Length, last, nameof(array));
            Data = array;
            Offset = offset;
            Rows = rows;
            Cols = cols;
            RowStride = rowStride;
            ColStride = colStride;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix(Vector vector)
        {
            Data = vector.Data;
            Offset = vector.Offset;
            Rows = vector.Length;
            RowStride = vector.Stride;
            Cols = 1;
            ColStride = Rows * vector.Stride;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix(Vector vector, nint rows, nint cols)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(vector.Length, rows * cols,
                nameof(vector));

            Data = vector.Data;
            Offset = vector.Offset;
            Rows = rows;
            RowStride = vector.Stride;
            Cols = cols;
            ColStride = rows * vector.Stride;
        }
        #endregion Constructors

        #region Alloc
        public static Matrix Create(nint rows, nint cols, bool uninited = false)
        {
            ArgumentOutOfRangeException.ThrowIfNegativeOrZero(rows, nameof(rows));
            ArgumentOutOfRangeException.ThrowIfNegativeOrZero(cols, nameof(cols));
            nint length = rows * cols;
            ArgumentOutOfRangeException.ThrowIfGreaterThan(length, Array.MaxLength, "size");
            double[] data = uninited ? GC.AllocateUninitializedArray<double>((int)length) 
                : new double[length];
            return new Matrix(data, 0, rows, cols);
        }

        public static Matrix Create(nint rows, nint cols, double val)
        {
            ArgumentOutOfRangeException.ThrowIfNegativeOrZero(rows, nameof(rows));
            ArgumentOutOfRangeException.ThrowIfNegativeOrZero(cols, nameof(cols));
            nint length = rows * cols;
            ArgumentOutOfRangeException.ThrowIfGreaterThan(length, Array.MaxLength, "size");
            double[] data = new double[length];
            Array.Fill(data, val);
            return new Matrix(data, 0, rows, cols, cols, 1);
        }
        #endregion Alloc

        #region Accessor and Slicer
        public ref double this[nint row, nint col]
        {
            get
            {
                if (row < 0 || row >= Rows)
                    throw new ArgumentOutOfRangeException(nameof(row), $"Row index {row} is out of range.");
                if (col < 0 || col >= Cols)
                    throw new ArgumentOutOfRangeException(nameof(col), $"Column index {col} is out of range.");
                return ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(Data), Offset + row * RowStride + col * ColStride);
            }
        }

        public Matrix this[NRange rows, NRange cols]
        {
            get
            {
                (nint startRow, nint lengthRow) = rows.GetOffsetAndLength(Rows);
                (nint startCol, nint lengthCol) = cols.GetOffsetAndLength(Cols);
                if (startRow < 0 || lengthRow < 0 || startRow + lengthRow > Rows)
                    throw new ArgumentOutOfRangeException(nameof(rows), $"Row range {rows} is out of range.");
                if (startCol < 0 || lengthCol < 0 || startCol + lengthCol > Cols)
                    throw new ArgumentOutOfRangeException(nameof(cols), $"Column range {cols} is out of range.");
                return SliceSubUncheck(startRow, lengthRow, startCol, lengthCol);
            }
            set
            {
                (nint startRow, nint lengthRow) = rows.GetOffsetAndLength(Rows);
                (nint startCol, nint lengthCol) = cols.GetOffsetAndLength(Cols);
                if (startRow < 0 || lengthRow < 0 || startRow + lengthRow > Rows)
                    throw new ArgumentOutOfRangeException(nameof(rows), $"Row range {rows} is out of range.");
                if (startCol < 0 || lengthCol < 0 || startCol + lengthCol > Cols)
                    throw new ArgumentOutOfRangeException(nameof(cols), $"Column range {cols} is out of range.");
                value.CopyTo(SliceSubUncheck(startRow, lengthRow, startCol, lengthCol));
            }
        }

        public Vector this[LengthIndex rowIndex, NRange cols]
        {
            get
            {
                nint row = rowIndex.GetOffset(Rows);
                (nint startCol, nint lengthCol) = cols.GetOffsetAndLength(Cols);
                if (row < 0 || row >= Rows)
                    throw new ArgumentOutOfRangeException(nameof(rowIndex),
                        $"Row index {row} is out of range.");
                if (startCol < 0 || lengthCol < 0 || startCol + lengthCol > Cols)
                    throw new ArgumentOutOfRangeException(nameof(cols),
                        $"Column range {cols} is out of range.");
                return SliceRowUncheck(row);
            }
            set
            {
                nint row = rowIndex.GetOffset(Rows);
                (nint startCol, nint lengthCol) = cols.GetOffsetAndLength(Cols);
                if (row < 0 || row >= Rows)
                    throw new ArgumentOutOfRangeException(nameof(rowIndex),
                        $"Row index {row} is out of range.");
                if (startCol < 0 || lengthCol < 0 || startCol + lengthCol > Cols)
                    throw new ArgumentOutOfRangeException(nameof(cols),
                        $"Column range {cols} is out of range.");
                value.CopyTo(SliceRowUncheck(row));
            }
        }

        public Vector this[NRange rows, LengthIndex colIndex]
        {
            get
            {
                (nint startRow, nint lengthRow) = rows.GetOffsetAndLength(Rows);
                nint col = colIndex.GetOffset(Cols);
                if (startRow < 0 || lengthRow < 0 || startRow + lengthRow > Rows)
                    throw new ArgumentOutOfRangeException(nameof(rows),
                        $"Row range {rows} is out of range.");
                if (col < 0 || col >= Cols)
                    throw new ArgumentOutOfRangeException(nameof(colIndex),
                        $"Column index {col} is out of range.");
                return SliceColUncheck(col, startRow, lengthRow);
            }
            set
            {
                (nint startRow, nint lengthRow) = rows.GetOffsetAndLength(Rows);
                nint col = colIndex.GetOffset(Cols);
                if (startRow < 0 || lengthRow < 0 || startRow + lengthRow > Rows)
                    throw new ArgumentOutOfRangeException(nameof(rows),
                        $"Row range {rows} is out of range.");
                if (col < 0 || col >= Cols)
                    throw new ArgumentOutOfRangeException(nameof(colIndex),
                        $"Column index {col} is out of range.");
                value.CopyTo(SliceColUncheck(col, startRow, lengthRow));
            }
        }

        public ref double GetHeadRef()
        {
            if(IsEmpty)
                throw new InvalidOperationException("Matrix is empty.");
            return ref Data[Offset];
        }

        public Vector GetRow(nint row)
        {
            if (row < 0 || row >= Rows)
                throw new ArgumentOutOfRangeException(nameof(row),
                    $"Row index {row} is out of range.");
            return SliceRowUncheck(row);
        }

        public Vector GetColumn(nint col)
        {
            if (col < 0 || col >= Cols)
                throw new ArgumentOutOfRangeException(nameof(col),
                    $"Column index {col} is out of range.");
            return SliceColUncheck(col);
        }

        public Vector Diag
        {
            get
            {
                nint diagLength = Math.Min(Rows, Cols);
                return new(Data, Offset, diagLength, RowStride + ColStride);
            }
            set
            {
                nint diagLength = Math.Min(Rows, Cols);
                var target = new Vector(Data, Offset, diagLength, RowStride + ColStride);
                BlasLike.Copy(value, target);
            }
        }

        internal ref double AtUncheck(nint row, nint col)
            => ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(Data),
                Offset + row * RowStride + col * ColStride);

        internal Vector SliceRowUncheck(nint row)
            => new(Data, (int)(Offset + row * RowStride), Cols, ColStride);

        internal Vector SliceRowUncheck(nint row, nint colStart)
            => new(Data, (int)(Offset + row * RowStride + colStart * ColStride),
                Cols - colStart, ColStride);

        internal Vector SliceRowUncheck(nint row, nint colStart, nint colLength)
            => new(Data, (int)(Offset + row * RowStride + colStart * ColStride),
                colLength, ColStride);

        internal Vector SliceColUncheck(nint col)
            => new(Data, (int)(Offset + col * ColStride), Rows, RowStride);

        internal Vector SliceColUncheck(nint col, nint rowStart)
            => new(Data, (int)(Offset + rowStart * RowStride + col * ColStride),
                Rows - rowStart, RowStride);

        internal Vector SliceColUncheck(nint col, nint rowStart, nint rowLength)
            => new(Data, (int)(Offset + rowStart * RowStride + col * ColStride),
                rowLength, RowStride);

        internal Matrix SliceSubUncheck(nint rowStart, nint colStart)
            => new(Data, (int)(Offset + rowStart * RowStride + colStart * ColStride),
                Rows - rowStart, Cols - colStart, RowStride, ColStride);

        internal Matrix SliceSubUncheck(nint rowStart, nint rowLength, nint colStart, nint colLength)
            => new(Data, (int)(Offset + rowStart * RowStride + colStart * ColStride),
                rowLength, colLength, RowStride, ColStride);
        #endregion Accessor and Slicer

        #region Memory
        public void FlattenTo(Span<double> target)
        {
            if (target.Length < TotalSize)
                throw new ArgumentException($"Target span is too small. " +
                    $"Required: {TotalSize}, Actual: {target.Length}");
            ref double targetRef = ref target[0];
            ref double head = ref Data[Offset];


            if (ColStride == 1)
            {
                for (int i = 0; i < Rows; i++)
                {
                    ExMemoryMarshal.CreateExSpan(ref head, (int)Cols)
                            .CopyTo(ExMemoryMarshal.CreateExSpan(ref targetRef, (int)Cols));
                    head = ref Unsafe.Add(ref head, RowStride);
                    targetRef = ref Unsafe.Add(ref targetRef, Cols);
                }
            }
            else
            {
                for (int i = 0; i < Rows; i++)
                {
                    ref double ptr = ref head;
                    int j = 0;
                    for (; j <= Cols - 4; j += 4)
                    {
                        targetRef = ptr;
                        Unsafe.Add(ref targetRef, 1) = Unsafe.Add(ref ptr, ColStride);
                        Unsafe.Add(ref targetRef, 2) = Unsafe.Add(ref ptr, ColStride * 2);
                        Unsafe.Add(ref targetRef, 3) = Unsafe.Add(ref ptr, ColStride * 3);
                        ptr = ref Unsafe.Add(ref ptr, ColStride * 4);
                        targetRef = ref Unsafe.Add(ref targetRef, 4);
                    }
                    for (; j < Cols; j++)
                    {
                        targetRef = Unsafe.Add(ref ptr, j * ColStride);
                        ptr = ref Unsafe.Add(ref ptr, ColStride);
                        targetRef = ref Unsafe.Add(ref targetRef, 1);
                    }
                    head = ref Unsafe.Add(ref head, RowStride);
                }
            }
        }

        public void Fill(double val)
        {
            BlasLike.Set(val, this);
        }

        public Matrix Transpose()
            => new(Data, Offset,
                rows: Cols, cols: Rows,
                rowStride: ColStride, colStride: RowStride);

        public void CopyTo(Matrix other)
        {
            BlasLike.Copy(this, other);
        }

        public Matrix Clone()
        {
            Matrix target = Create(Rows, Cols, true);
            BlasLike.Copy(this, target);
            return target;
        }

        public Matrix MakeContinous(bool forceCopy)
        {
            if (!forceCopy && RowStride == Cols && ColStride == 1)
            {
                return this; // Already continuous
            }
            return Clone();
        }

        public Vector Flatten()
        {
            double[] arr = new double[Length];
            FlattenTo(arr);
            return new(arr);
        }

        internal Span<double> GetSpan()
        {
            if (IsEmpty)
                throw new InvalidOperationException("Matrix is empty.");
            return Data.AsSpan(Offset);
        }

        /// <summary>
        /// If the matrix is not row-major layer,
        /// it will create a new continuous matrix.
        /// </summary>
        /// <remarks>Just use for BLAS.</remarks>
        /// <returns>A row-majo</returns>
        public Matrix MakeRowMajor()
        {
            if (ColStride == 1)
            {
                return this; // Already continuous
            }
            return Clone();
        }
        #endregion Memory

        #region IEnumerable
        IEnumerator<double> IEnumerable<double>.GetEnumerator()
            => new MatrixEnumberator(this);

        IEnumerator IEnumerable.GetEnumerator()
            => new MatrixEnumberator(this);

        public RefMatrixEnumberator GetEnumerator()
            => new(this);

        internal struct MatrixEnumberator : IEnumerator<double>
        {
            private Matrix _matrix;
            private nint _currentRow;
            private nint _currentCol;
            public MatrixEnumberator(Matrix matrix)
            {
                _matrix = matrix;
                _currentRow = 0;
                _currentCol = 0;
            }

            public readonly double Current
                => _matrix.AtUncheck(_currentRow, _currentCol);

            readonly double IEnumerator<double>.Current => Current;

            readonly object IEnumerator.Current => Current;

            public bool MoveNext()
            {
                if(_currentCol < _matrix.Cols - 1)
                {
                    _currentCol++;
                    return true;
                }
                else if (_currentRow < _matrix.Rows - 1)
                {
                    _currentCol = 0;
                    _currentRow++;
                    return true;
                }
                else
                {
                    return false;
                }
            }

            public void Reset()
            {
                _currentRow = 0;
                _currentCol = 0;
            }

            void IDisposable.Dispose() { }
        }

        public ref struct RefMatrixEnumberator
        {
            private ref double _currentVal;
            private nint _currentRow;
            private readonly nint _rowEnd;
            private readonly nint _rowStride;
            private nint _currentCol;
            private readonly nint _colEnd;
            private readonly nint _colStride;
            private bool _init;
            public RefMatrixEnumberator(Matrix matrix)
            {
                _currentRow = 0;
                _currentCol = 0;
                _rowEnd = matrix.Rows - 1;
                _colEnd = matrix.Cols - 1;
                _rowStride = matrix.RowStride;
                _colStride = matrix.ColStride;
                _currentVal = ref matrix.AtUncheck(0, 0);
                _init = false;
            }
            public ref double Current
                => ref _currentVal;

            public bool MoveNext()
            {
                if (!_init)
                {
                    _init = true;
                    return _currentRow <= _rowEnd && _currentCol <= _colEnd;
                }

                if (_currentCol < _colEnd)
                {
                    _currentCol++;
                    _currentVal = ref Unsafe.Add(ref _currentVal, _colStride);
                    return true;
                }
                else if (_currentRow < _rowEnd)
                {
                    _currentCol = 0;
                    _currentVal = ref Unsafe.Subtract(ref _currentVal, _colEnd * _colStride);
                    _currentRow++;
                    _currentVal = ref Unsafe.Add(ref _currentVal, _rowStride);
                    return true;
                }
                else
                    return false;

            }

            public void Reset()
            {
                _currentVal = ref Unsafe.Subtract(ref _currentVal,
                    _currentRow * _colStride + _currentCol * _rowStride);
                _currentRow = 0;
                _currentCol = 0;
            }
        }

        #endregion IEnumerable

        #region Math
        
        public static Matrix operator +(Matrix left, Matrix right)
        {
            right = right.Clone();
            BlasLike.Add(left, right);
            return right;
        }

        public static Matrix operator +(double left, Matrix right)
        {
            right = right.Clone();
            var rightView = right.Clone().AsSNTensor();

            System.Numerics.Tensors.Tensor.Add(rightView, left);
            return right;
        }

        public static Matrix operator +(Matrix left, double right)
        {
            left = left.Clone();
            var leftView = left.Clone().AsSNTensor();

            System.Numerics.Tensors.Tensor.Add(leftView, right);
            return left;
        }

        public static Matrix operator -(Matrix left, Matrix right)
        {
            left = left.Clone();
            BlasLike.Sub(right, left);
            return left;
        }

        public static Matrix operator -(Matrix left, double right)
        {
            left = left.Clone();
            var leftView = left.Clone().AsSNTensor();

            System.Numerics.Tensors.Tensor.Subtract(leftView, right);
            return left;
        }

        public static Matrix operator -(double left, Matrix right)
        {
            right = right.Clone();
            var rightView = right.Clone().AsSNTensor();

            System.Numerics.Tensors.Tensor.Subtract(rightView, left);
            return right;
        }

        public static Matrix operator *(double left, Matrix right)
        {
            var dest = Create(right.Rows, right.Cols, uninited:true);
            BlasLike.Scal2(left, right, dest);
            return dest;
        }

        public static Matrix operator *(Matrix left, double right)
            => right * left;

        public static Matrix operator /(Matrix left, double right)
            => left * (1.0 / right);

        public static Vector operator *(Matrix left, Vector right)
        {
            var dest = Vector.Create(left.Rows, true);
            BlasLike.Gemv(1.0, left, right, 0.0, dest);
            return dest;
        }

        public static Matrix operator *(Matrix left, Matrix right)
        {
            ArgumentOutOfRangeException.ThrowIfNotEqual(left.Cols, right.Rows,
                "Matrix multiplication requires the number of columns in the left matrix " +
                "to be equal to the number of rows in the right matrix.");
            nint m = left.Rows;
            nint n = right.Cols;
            nint k = left.Cols;
            var dest = Matrix.Create(m, n, uninited: true);
            left = left.MakeRowMajor();
            right = right.MakeRowMajor();
            checked
            {
                Blas.gemm(Layout.RowMajor,
                    TransA: Trans.No, TransB: Trans.No,
                    M: (int)m, N: (int)n, K: (int)k,
                    alpha: 1.0,
                    left.GetSpan(), (int)left.RowStride,
                    right.GetSpan(), (int)right.RowStride,
                    beta: 0.0,
                    dest.GetSpan(), (int)m
                    );
            }
            return dest;
        }

        #endregion math

        #region lapack

        public SEvd SEvd() => new(this);

        public LU LU(bool inplace = false) => new(this, inplace);

        #endregion lapack

        #region ufunc

        public Matrix Pointwise<TUnaryAction>()
            where TUnaryAction : struct, IUnaryOperator<double, double>
        {
            Matrix dest = Clone();
            BlasLike.ApplyM<TUnaryAction>(dest);
            return dest;
        }

        public Matrix PointwiseAbs()
            => Pointwise<AbsOperator<double>>();


        #endregion ufunc

        public struct DMatrixSpanDebugView
        {
            public object Items;
            public DMatrixSpanDebugView(Matrix span)
            {
                Rows = span.Rows;
                Cols = span.Cols;
                RowStride = span.RowStride;
                ColStride = span.ColStride;

                if (span.IsEmpty || TooLong)
                {
                    Items = Array.Empty<double>();
                }
                else
                {
                    Dictionary<(long x, long y), double>
                        values = new((int)span.TotalSize);
                    for (nint i = 0; i < span.Rows; i++)
                    {
                        for (nint j = 0; j < span.Cols; j++)
                        {
                            values[(i, j)] = span.AtUncheck(i, j);
                        }
                    }
                    Items = values;
                }
            }

            public long Rows { get; }
            public long Cols { get; }
            public long RowStride { get; }
            public long ColStride { get; }

            [DebuggerBrowsable(DebuggerBrowsableState.Never)]
            public bool TooLong => Rows * Cols > TooLongLimit;

            [DebuggerBrowsable(DebuggerBrowsableState.Never)]
            public const int TooLongLimit = int.MaxValue / 2;
        }
    }
}
