using MKLNET;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SimpleHelpers.LinearAlg.Dec
{
    public class SEvd
    {
        readonly Vector eigenValues;
        readonly Matrix eigenVectors;
        readonly Matrix origin;
        readonly bool lower;
        readonly nint n;
        bool computed;
        UpLoChar Uplo => lower ? UpLoChar.Lower : UpLoChar.Upper;

        public Vector EigenValues
        {
            get
            {
                Compute();
                return eigenValues;
            }
        }

        public Matrix EigenVectors
        {
            get
            {
                Compute();
                return eigenVectors;
            }
        }

        public SEvd(Matrix mat, bool lower = true, bool inplace = false)
        {
            if (mat.Rows > int.MaxValue)
                throw new NotImplementedException("Matrix size exceeds int.MaxValue, not supported yet.");
            ArgumentOutOfRangeException.ThrowIfNotEqual(mat.Rows, mat.Cols, nameof(mat));
            n = mat.Rows;
            origin = mat;
            this.lower = lower;
            eigenValues = Vector.Create(n, uninited: true);
            if (inplace)
            {
                if (mat.ColStride != 1)
                    throw new ArgumentException("In-place operation " +
                        "requires row-major matrix.", nameof(mat));
                eigenVectors = mat;
            }
            else
            {
                eigenVectors = mat.Clone();
            }
            computed = false;
        }

        private void Compute()
        {
            if (computed)
                return;
            int info = Lapack.syev(
                Layout.RowMajor,
                'V', uplo:Uplo,
                (int)n,
                eigenVectors.GetSpan(), (int)n, eigenValues.GetSpan());
            if (info != 0)
                throw new LinalgException("dsyev", info,
                    $"LAPACK dsyev failed.");
            computed = true;
        }

        public Matrix Power(double n)
        {
            var diag = EigenValues.Clone();
            BlasLike.Pow(n, diag);
            return ReconstructMatrix(diag);
        }

        public Matrix Sqrt()
        {
            var diag = EigenValues.Clone();
            BlasLike.Sqrt(diag);
            return ReconstructMatrix(diag);
        }

        public Matrix ReverseSqrt()
        {
            var diag = EigenValues.Clone();
            for (int i = 0; i < n; i++)
            {
                if (diag[i] <= 0)
                    throw new LinalgException("SEVD.ReverseSqrt",
                        $"Eigenvalue {i} is non-positive: {diag[i]}.");
                diag[i] = 1.0 / Math.Sqrt(diag[i]);
            }
            return ReconstructMatrix(diag);
        }

        public Matrix ReconstructMatrix(Vector diag)
        {
            Matrix diagMat = Matrix.Create(n, n);
            diagMat.Diag = diag;
            return eigenVectors * diagMat * eigenVectors.Transpose();
        }
    }
}
