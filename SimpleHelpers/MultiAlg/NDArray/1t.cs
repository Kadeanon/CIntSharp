using NetFabric.Numerics.Tensors.Operators;
using SimpleHelpers.LinearAlg;
using Details = SimpleHelpers.LinearAlg.BlasLike.Details;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {
        public static void Add(NDArray x, NDArray y)
            => ApplyTo<AddOperator<double>>(x, y);

        public static void Axpy(double alpha, NDArray x, NDArray y)
            => ApplyToWith<Details.AxpyOperator, double>(x, alpha, y);

        public static void Axpy(Vector alpha, NDArray x, NDArray y)
            => BatchApplyToWith<Details.AxpyOperator>(x, alpha, y);

        public static void Copy(NDArray x, NDArray y)
            => ApplyAssign<IdentityOperator<double>>(x, y);

        public static void Invert(NDArray x)
            => ApplyWith<Details.InvertOperator, double>(x, 1);

        public static void InvScal(double alpha, NDArray x)
        {
            if (alpha == 0.0 || alpha == -0.0 || !double.IsFinite(alpha))
            {
                throw new ArgumentException("Error: alpha must be a finite non-zero value.");
            }
            ApplyWith<MultiplyOperator<double>, double>(x, 1 / alpha);
        }

        public static void Scal(double alpha, NDArray x)
            => ApplyWith<MultiplyOperator<double>, double>(x, alpha);

        public static void Scal2(double alpha, NDArray x, NDArray y)
            => ApplyWithAssign<MultiplyOperator<double>, double>(x, alpha, y);

        public static void Pow(double alpha, NDArray x)
            => ApplyWith<Details.PowOperator, double>(x, alpha);

        public static void Sqrt(NDArray x)
            => Apply<SqrtOperator<double>>(x);

        public static void Set(double alpha, NDArray x)
            => ApplySelf<IdentityOperator<double>, double>(x, alpha);

        public static void BatchSet(Vector alpha, NDArray y)
            => BatchApplySelf<IdentityOperator<double>>(y, alpha);

        public static void Sub(NDArray x, NDArray y)
            => ApplyTo<Details.SubtractedByOperator>(x, y);

        public static void Xpby(NDArray x, double beta, NDArray y)
            => ApplyToWith<Details.XpbyOperator, double>(x, beta, y);

        public static void Xpby(NDArray x, Vector beta, NDArray y)
            => BatchApplyToWith<Details.XpbyOperator>(x, beta, y);
    }
}
