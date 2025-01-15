using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Native
{
    [AttributeUsage(AttributeTargets.Assembly)]
    public class LibraryPathAttribute(string path) : Attribute
    {
        public string Path { get; } = path;

        public static LibraryPathAttribute? Instance =>
            Assembly.GetExecutingAssembly().GetCustomAttribute<LibraryPathAttribute>();

        public string GetDllPlatformPath(string name) 
        {
            string path = System.IO.Path.Join(this.Path, "Native", name);
            var platform = GetPlatform();
            if (platform == OSPlatform.Windows)
            {
                path = System.IO.Path.Combine(path, "Windows-x64", $"{name.ToLower()}.dll");
            }
            return path;
        }

        public static OSPlatform GetPlatform()
        {
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                return OSPlatform.Windows;
            /*
            else if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                platformKey = OSPlatform.Linux;
            else if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                platformKey = OSPlatform.OSX;
            */
            else
                throw new PlatformNotSupportedException("不支持的操作系统平台。");
        }
    }
}
