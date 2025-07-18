using MKLNET;

namespace SimpleHelpers.LinearAlg.Decs
{
    public class LU
    {
        readonly Matrix matrix;
        readonly int[] ipiv;
        readonly int n;
        bool computed;

        public LU(Matrix mat, bool inplace = false)
        {
            if (mat.Rows > int.MaxValue)
                throw new NotImplementedException("Matrix size exceeds int.MaxValue, not supported yet.");
            ArgumentOutOfRangeException.ThrowIfNotEqual(mat.Rows, mat.Cols, nameof(mat));
            n = (int)mat.Rows;
            if (inplace)
            {
                if (mat.ColStride != 1)
                    throw new ArgumentException("In-place operation " +
                        "requires row-major matrix.", nameof(mat));
                matrix = mat;
            }
            else
            {
                matrix = mat.Clone();
            }
            ipiv = new int[n];
            computed = false;
        }

        private void Compute()
        {
            if (computed)
                return;
            var info = Lapack.getrf(Layout.RowMajor,
                n, n, matrix.GetSpan(), n, ipiv);
            if (info != 0)
            {
                if (info > 0)
                {
                    int index = info - 1;
                    throw new LinalgException("dgetrf", info,
                        "LU decomposition failed, " +
                        $"the {index} diagonal element of matrix is zero, " +
                        $"and the solve could not be completed.");
                }
                else
                {
                    throw new LinalgException("dgetrf", info,
                        $"LAPACK getrf failed.");
                }
            }
            computed = true;
        }

        public Vector Solve(Vector b)
        {
            Compute();
            if (b.Length != n)
                throw new ArgumentException("RHS vector size mismatch", nameof(b));
            var x = b.Clone();
            var info = Lapack.getrs(Layout.RowMajor, TransChar.No,
                n, nrhs: 1,
                matrix.GetSpan(), (int)matrix.RowStride,
                ipiv, 
                x.GetSpan(), 1);
            if (info != 0)
                throw new LinalgException("dgetrs", info, 
                    $"LAPACK getrs failed with info = {info}");

            return x;
        }

        public Matrix Solve(Matrix b)
        {
            Compute();
            if (b.Rows != n)
                throw new ArgumentException("RHS matrix row size mismatch", nameof(b));
            int nrhs = (int)b.Cols;
            b = b.MakeRowMajor();
            var info = Lapack.getrs(Layout.RowMajor, TransChar.No,
                n, nrhs,
                matrix.GetSpan(), (int)matrix.RowStride,
                ipiv, 
                b.GetSpan(), (int)b.RowStride);
            if (info != 0)
                throw new LinalgException("dgetrs", info,
                    $"LAPACK getrs failed.");
            return b;
        }
    }
}
