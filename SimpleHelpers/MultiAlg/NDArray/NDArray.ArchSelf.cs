using NetFabric.Numerics.Tensors.Operators;
using SimpleHelpers.LinearAlg;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {

        public NDArray ScaledBy(double alpha)
        {
            if(alpha == 0.0)
                Fill(0.0);
            else 
            if (alpha != 1.0)
                ApplyWith<MultiplyOperator<double>, double>(this, alpha);
            return this;
        }

        public NDArray AddedBy(double alpha)
        {
            if (alpha != 0.0)
                ApplyWith<AddOperator<double>, double>(this, alpha);
            return this;
        }

        public NDArray AddedBy(NDArray other)
        {
            ArgumentNullException.ThrowIfNull(other, nameof(other));
            ApplyTo<AddOperator<double>>(other, this);
            return this;
        }

        public NDArray AddedByScaled(double factor, NDArray other)
        {
            ArgumentNullException.ThrowIfNull(other, nameof(other));
            ApplyToWith<BlasLike.Details.AxpyOperator, double>(other, factor, this);
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
            if (alpha != 0.0)
                ApplyWith<SubtractOperator<double>, double>(this, alpha);
            return this;
        }

        public NDArray AssignedBy(NDArray other)
        {
            ApplyAssign<IdentityOperator<double>>(other, this);
            return this;
        }
    }
}
