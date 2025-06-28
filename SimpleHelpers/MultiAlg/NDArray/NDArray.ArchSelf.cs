using NetFabric.Numerics.Tensors.Operators;
using SimpleHelpers.LinearAlg;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {

        public NDArray ScaledBy(double alpha)
        {
            if (alpha == 1.0)
            {
                return this;
            }
            ApplyWith<MultiplyOperator<double>, double>(this, alpha);
            return this;
        }

        public NDArray AddedBy(double alpha)
        {
            if (alpha == 0.0)
            {
                return this;
            }
            ApplyWith<AddOperator<double>, double>(this, alpha);
            return this;
        }

        public NDArray AddedBy(NDArray other)
        {
            ArgumentNullException.ThrowIfNull(other, nameof(other));
            ApplyTo<AddOperator<double>>(other, this);
            return this;
        }

        public NDArray SubtractedBy(NDArray other)
        {
            ArgumentNullException.ThrowIfNull(other, nameof(other));
            ApplyTo<BlasLike.Details.SubtractedByOperator>(other, this);
            return this;
        }

        public NDArray SubtractedBy(double alpha)
        {
            if (alpha == 0.0)
            {
                return this;
            }
            ApplyWith<SubtractOperator<double>, double>(this, alpha);
            return this;
        }
    }
}
