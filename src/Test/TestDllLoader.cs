using CintSharp.Native.Libcint;

namespace CintSharp.Test
{
    [TestClass]
    public sealed class TestDllLoader
    {
        [TestMethod]
        public void TestLibcint()
        {
            LibcintHandler.LoadDll();
            Assert.IsNotNull(LibcintHandler.Instance);
            Assert.IsFalse(LibcintHandler.Instance.IsInvalid);
        }


        [TestMethod]
        public void TestInvalidLibcint()
        {
            Assert.ThrowsException<DllNotFoundException>(() => 
                LibcintHandler.LoadDll("libcint.dll"));
        }
    }
}
