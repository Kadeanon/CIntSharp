using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SimpleHelpers.LinearAlg
{
    internal class LinalgException : Exception
    {
        internal string FunctionName { get; }

        internal int ErrorCode { get; }

        public LinalgException(string functionName, 
            int errorCode, string message) : base(
                $"MKL function {functionName} failed " +
            $"with errorCode {errorCode}: {message}")
        {
            FunctionName = functionName;
            ErrorCode = errorCode;
        }

        public LinalgException(string functionName, string message)
            : base($"Managed Linalg function {functionName} " +
                  $"failed: {message}")
        {
            FunctionName = functionName;
            ErrorCode = 0;
        }
    }
}
