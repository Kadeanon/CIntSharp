using SimpleHelpers;
using SimpleHelpers.LinearAlg;
using SimpleHelpers.MultiAlg;

namespace Examples
{
    public abstract class DIISBase
        (nint errorLength, int maxDiis)
    {
        protected readonly int maxDiis = maxDiis;
        protected readonly Matrix errorOvlp = 
            Matrix.Create(maxDiis + 1, maxDiis + 1);
        protected readonly Matrix errors = 
            Matrix.Create(maxDiis, errorLength);
        protected int currentSize = 0;
        protected int currentIndex = -1;

        protected Vector UpdateAndGetCoeffs(Vector error)
        {
            currentIndex = (currentIndex + 1) % maxDiis;
            currentSize = Math.Clamp(currentSize + 1, 0, maxDiis);
            this.errors[currentIndex, ..] = error;

            if (currentSize <= 1)
                return Vector.Create(currentSize);

            var workSize = currentSize + 1;
            var errors = errorOvlp[..workSize, ..workSize];
            for (int i = 0; i < currentSize; i++)
            {
                var error_i = this.errors[i, ..];
                var ovlp = error * error_i;
                errors[i, currentSize] = -1;
                errors[i, currentIndex] = ovlp;
                errors[currentIndex, i] = ovlp;
                errors[currentSize, i] = -1;
            }
            errors[currentSize, currentSize] = 0;
            var b = Vector.Create(workSize);
            b[currentSize] = -1;
            var lu = errors.LU();

            var coeffs = lu.Solve(b);
            coeffs[currentIndex] -= 1;

            return coeffs[..currentSize];
        }

        protected Vector UpdateAndGetCoeffs(Action<Vector> setter)
        {
            currentIndex = (currentIndex + 1) % maxDiis;
            currentSize = Math.Clamp(currentSize + 1, 0, maxDiis);
            var error = this.errors.GetRow(currentIndex);
            setter(error);

            if (currentSize <= 1)
                return Vector.Create(currentSize);

            var workSize = currentSize + 1;
            var errors = errorOvlp[..workSize, ..workSize];
            for (int i = 0; i < currentSize; i++)
            {
                var error_i = this.errors[i, ..];
                var ovlp = error * error_i;
                errors[i, currentSize] = -1;
                errors[i, currentIndex] = ovlp;
                errors[currentIndex, i] = ovlp;
                errors[currentSize, i] = -1;
            }
            errors[currentSize, currentSize] = 0;
            var b = Vector.Create(workSize);
            b[currentSize] = -1;
            var lu = errors.LU();

            var coeffs = lu.Solve(b);
            coeffs[currentIndex] -= 1;

            return coeffs[..currentSize];
        }
    }
}
