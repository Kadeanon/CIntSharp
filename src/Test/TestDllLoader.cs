using CintSharp.Native;
using CintSharp.Native.Libcint;

namespace CintSharp.Test
{
    [TestClass]
    public sealed class TestDllLoader
    {
        [TestMethod]
        public void TestLibcint()
        {
            Assert.IsNotNull(LibcintHandler.Instance);
            Assert.IsFalse(LibcintHandler.Instance.IsInvalid);
        }
    }
}
