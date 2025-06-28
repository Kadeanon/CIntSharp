using NetFabric.Numerics.Tensors.Operators;
using System.Numerics.Tensors;

namespace SimpleHelpers.MultiAlg
{
    public partial class NDArray
    {
        public static NDArray operator +(NDArray left, NDArray right)
        {
            var result = right.Clone();
            ApplyTo<AddOperator<double>>(left, result);
            return result;
        }

        public static NDArray operator -(NDArray left, NDArray right)
        {
            var result = right.Clone();
            ApplyTo<SubtractOperator<double>>(left, result);
            return result;
        }

        public static NDArray operator +(NDArray left, double right)
        {
            var result = left.Clone();
            ApplyWith<AddOperator<double>, double>(result, right);
            return result;
        }

        public static NDArray operator +(double left, NDArray right)
            => right + left;

        public static NDArray operator -(NDArray left, double right)
        {
            var result = left.Clone();
            ApplyWith<SubtractOperator<double>, double>(result, right);
            return result;
        }

        public static NDArray operator -(double left, NDArray right)
        {
            throw new NotImplementedException("Addition of two NDArray instances is not implemented yet.");
        }

        public static NDArray operator *(NDArray left, double right)
        {
            var result = left.Clone();
            ApplyWith<MultiplyOperator<double>, double>(result, right);
            return result;
        }

        public static NDArray operator *(double left, NDArray right)
            => right * left;

        public static NDArray operator /(NDArray left, double right)
        {
            var result = left.Clone();
            ApplyWith<DivideOperator<double>, double>(result, right);
            return result;
        }

        public double Norm()
        {
             return Tensor.Norm(Clone().AsReadOnlyTensorSpan());
        }
    }
}
