using Microsoft.Win32.SafeHandles;
using System;
using System.Collections.Concurrent;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Native
{
    public class DllHandler : SafeHandleZeroOrMinusOneIsInvalid
    {
        /// <summary>
        /// Dictionary to store the libcint dll. Default is empty string to local default direction.
        /// </summary>
        public static string LibcintPath { get; set; } = "";
        public ConcurrentDictionary<string, Delegate> ApiDict { get; }

        public DllHandler(string dllName) : base(true)
        {
            Console.WriteLine($"Try load dll: {dllName}......");
            try
            {
                handle = NativeLibrary.Load(dllName,
                    Assembly.GetExecutingAssembly(),
                    DllImportSearchPath.LegacyBehavior);
            }
            catch (DllNotFoundException ex) 
            {
                Console.WriteLine($"Fail load - dll not found: {ex.Message}");
                throw;
            }
            catch (BadImageFormatException ex) 
            {
                Console.WriteLine($"Fail load - bad image format: {ex.Message}");
                throw;
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Fail load - unknown error: {ex.Message}");
                throw;
            }
            Console.WriteLine($"Load successfully!");
            ApiDict = [];
        }

        protected override bool ReleaseHandle()
        {
            ApiDict.Clear();
            NativeLibrary.Free(handle);
            handle = IntPtr.Zero;
            return true;
        }

        protected override void Dispose(bool disposing)
        {
            if (!this.IsInvalid)
            {
                ReleaseHandle();
            }
        }

        /// <summary>
        /// Convert the function pointer to a delegate and store it in the dictionary.
        /// </summary>
        /// <typeparam name="TDelegate"> the type of the delegate to convert to</typeparam>
        /// <param name="apiName"> the name of function </param>
        /// <returns> the delegate converted from the function </returns>
        /// <exception cref="ArgumentException"></exception>
        /// <exception cref="System.ComponentModel.Win32Exception"></exception>
        public TDelegate Invoke<TDelegate>(string apiName) where TDelegate : Delegate
        {
            var api = ApiDict.GetOrAdd(apiName,
                name =>
                {
                    if (!NativeLibrary.TryGetExport(handle, name, out var apiPointer))
                    {
                        throw new System.ComponentModel.Win32Exception(Marshal.GetLastWin32Error(), name);
                    }
                    return Marshal.GetDelegateForFunctionPointer<TDelegate>(apiPointer);
                });
            return api as TDelegate ??
                throw new ArgumentException("The ApiName is already in the dictionary, " +
                    "but the type is not the same as the type of the delegate to convert to." +
                    $"The api is a {api.GetType().Name} but is expected as a {typeof(TDelegate).Name} delegate.");
        }
    }
}
