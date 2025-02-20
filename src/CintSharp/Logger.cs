using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp
{
    internal static class Logger
    {

        public static void Log(string message, params object[] objects)
        {
            Console.WriteLine(message, objects);
        }

        public static void Error(string message, params object[] objects)
        {
            Console.Write("Error:");
            Console.WriteLine(message, objects);
        }

    }
}
