using System.Buffers;

namespace SimpleHelpers.MultiAlg.TensorContract
{
    public static partial class ContractMethods
    {
        public static ArrayPool<nint> ContractNintPool { get; } = ArrayPool<nint>.Create();
    }
}
