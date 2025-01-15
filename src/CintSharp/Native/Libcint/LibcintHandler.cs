using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Native.Libcint
{
    public class LibcintHandler : DllHandler
    {
        public LibcintHandler(string path) : base(path)
        {
        }

        public static string Name => "Libcint";

        public static LibcintHandler? Instance { get; private set; }

        internal const string DefaultDllPath = @"/Native/Libcint/Windows-x64/libcint.dll";

        public static void LoadDll(string path = "") 
        {
            if(String.IsNullOrWhiteSpace(path))
            {
                var attr = LibraryPathAttribute.Instance ?? throw new InvalidOperationException("The LibraryPathAttribute is not found in the assembly.");
                path = attr.GetDllPlatformPath(Name);
            }
            Instance?.Dispose();
            Instance = new LibcintHandler(path);
        }
    }
}
