using Microsoft.Win32.SafeHandles;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.Native
{
    public class DllHandler : SafeHandleZeroOrMinusOneIsInvalid
    {

        public Dictionary<string, Delegate> ApiDict { get; }

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
        /// <param name="ApiName"> the name of function </param>
        /// <returns> the delegate converted from the function </returns>
        /// <exception cref="ArgumentException"></exception>
        /// <exception cref="System.ComponentModel.Win32Exception"></exception>
        public TDelegate Invoke<TDelegate>(string ApiName) where TDelegate : Delegate
        {
            if (ApiDict.TryGetValue(ApiName, out var api))
            {
                if (api is TDelegate dele)
                {
                    return dele;
                }
                else
                {
                    throw new ArgumentException("The ApiName is already in the dictionary, " +
                        "but the type is not the same as the type of the delegate to convert to." +
                        $"The api is a {api.GetType().Name} but is expected as a {typeof(TDelegate).Name} delegate.");
                }
            }
            else
            {
                if (!NativeLibrary.TryGetExport(handle, ApiName, out var apiPointer))
                {
                    throw new System.ComponentModel.Win32Exception(Marshal.GetLastWin32Error(), ApiName);
                }
                var apiDele = Marshal.GetDelegateForFunctionPointer<TDelegate>(apiPointer);
                ApiDict.Add(ApiName, apiDele);
                return apiDele;
            }
        }
    }
}
