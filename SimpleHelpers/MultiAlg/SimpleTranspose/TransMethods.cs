using NetFabric.Numerics.Tensors.Operators;

namespace SimpleHelpers.MultiAlg.SimpleTranspose
{
    public class TransMethods
    {
        public static void CopyTo(NDArray source, NDArray destination)
            => NDArray.ApplyAssign<IdentityOperator<double>>(source, destination);
    }
}
