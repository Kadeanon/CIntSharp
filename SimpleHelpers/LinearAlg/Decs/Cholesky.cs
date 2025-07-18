using MKLNET;
using SimpleHelpers.MultiAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SimpleHelpers.LinearAlg.Decs
{
    public class Cholesky
    {
        readonly Matrix matrix;
        readonly int n;
        readonly bool lower;
        bool computed;
        UpLoChar Uplo => lower ? UpLoChar.Lower : UpLoChar.Upper;


        public Cholesky(Matrix mat, bool inplace = false, bool lower = false)
        {
            if (mat.Rows > int.MaxValue)
                throw new NotImplementedException("Matrix size exceeds int.MaxValue, not supported yet.");
            ArgumentOutOfRangeException.ThrowIfNotEqual(mat.Rows, mat.Cols, nameof(mat));
            n = (int)mat.Rows;
            if(inplace)
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
            computed = false;
            this.lower = lower;
        }

        private void Compute()
        {
            if (computed)
                return;
            var info = Lapack.potrf(Layout.RowMajor,
                Uplo, n, matrix.GetSpan(), (int)matrix.RowStride);
            if (info != 0)
            {
                if (info > 0)
                {
                    int index = info - 1;
                    throw new LinalgException("dpotrf", info,
                        "Chol decomposition failed, the leading minor " +
                        $"of order {index} is not positive-definite, " +
                        "and the solve could not be completed.");
                }
                else
                {
                    throw new LinalgException("dpotrf", info,
                        $"LAPACK dpotrf failed.");
                }
            }
            computed = true;
        }

        public Vector Solve(Vector b, bool inplace = false)
        {
            Compute();
            if (b.Length != n)
                throw new ArgumentException("RHS vector size mismatch", nameof(b));
            var x = b.Clone();
            var info = Lapack.potrs(Layout.RowMajor, Uplo,
                n, nrhs: 1, 
                matrix.GetSpan(), (int)matrix.RowStride, 
                x.GetSpan(), 1);
            if (info != 0)
                throw new LinalgException("dpotrs", info,
                    $"LAPACK potrs failed with info = {info}");

            return x;
        }

        public Matrix Solve(Matrix b, bool inplace = false)
        {
            Compute();
            if (b.Rows != n)
                throw new ArgumentException("RHS matrix row size mismatch", nameof(b));
            int nrhs = (int)b.Cols;
            b = b.Clone();
            var info = Lapack.potrs(Layout.RowMajor, Uplo,
                n, nrhs, 
                matrix.GetSpan(), (int)matrix.RowStride, 
                b.GetSpan(), (int)b.RowStride);
            if (info != 0)
                throw new LinalgException("dpotrs", info,
                    $"LAPACK potrs failed.");
            return b;
        }
    }
}
