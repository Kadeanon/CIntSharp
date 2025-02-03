using CintSharp.DataStructures;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.BasisParser
{
    public interface IBasisParser
    {
        Shells ParseBasis(BasisDict.BasisKey key, List<double> envs);
        bool TryParseBasis(BasisDict.BasisKey key, List<double> envs, [NotNullWhen(true)]out Shells? entry) 
        {
            try 
            {
                entry = ParseBasis(key, envs);
                return true;
            }
            catch
            {
                entry = null;
                return false;
            }
        }
    }
}
