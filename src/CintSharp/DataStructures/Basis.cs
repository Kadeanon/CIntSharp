global using Shells = System.Collections.Generic.List<CintSharp.DataStructures.Shell>;
using CintSharp.BasisParser;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO.Pipes;
using System.Linq;
using System.Text;
using System.Text.Json.Nodes;
using System.Threading.Tasks;


namespace CintSharp.DataStructures
{
    /// <summary>
    /// A dictionary of basis sets.
    /// </summary>
    /// <param name="parser">A <see cref="IBasisParser"/> instance to parse basis sets.</param>
    /// <param name="rawList"></param>
    public class BasisDict(IBasisParser parser, List<double> rawList)
    {

        internal IBasisParser DataParser { get; } = parser;

        internal List<double> RawList { get; } = rawList;

        internal Dictionary<BasisKey, Shells> RawDict { get; } = [];

        public Shells this[BasisKey key]
        {
            get => Get(key);
            set => Set(key, value);
        }

        public Shells this[string basisName, int atomNumber]
        {
            get => Get(new(basisName, atomNumber));
            set => Set(new(basisName, atomNumber), value);
        }

        public bool TryAdd(BasisKey key, [NotNullWhen(true)] out Shells? entry) 
        {
            if (RawDict.TryGetValue(key, out entry))
            {
                return true;
            }

            try
            {
                entry = AddUnchecked(key);
                return true;
            }
            catch (Exception) 
            {
                entry = null;
                return false; 
            }
        }

        public Shells Add(BasisKey key)
        {
            if (RawDict.TryGetValue(key, out Shells? entry))
            {
                return entry;
            }
            try
            {
                return AddUnchecked(key);
            }
            catch (Exception e)
            {
                throw new ArgumentException
                    ($"Failed to add basis {key.BasisName} for atom {key.AtomNumber}", e);
            }
        }

        public bool ContainsKey(BasisKey key)
        {
            return RawDict.ContainsKey(key);
        }

        public bool ContainsKey(string basisName, int atomNumber)
        {
            return RawDict.ContainsKey(new(basisName, atomNumber));
        }

        Shells Get(BasisKey key)
        {
             if (RawDict.TryGetValue(key, out Shells? entry))
            {
                return entry;
            }
            try 
            {
                entry = DataParser.ParseBasis(key, RawList);
                RawDict.Add(key, entry);
                return entry;
            }
            catch(Exception e) 
            {
                throw new ArgumentException($"Failed to get basis {key.BasisName} for atom {key.AtomNumber}", e);
            }
        }

        void Set(BasisKey key, Shells entry)
        {
            RawDict[key] = entry;
        }

        Shells AddUnchecked(BasisKey key)
        {
            var entry = DataParser.ParseBasis(key, RawList);
            RawDict.Add(key, entry);
            return entry;
        }

    }
    public record struct BasisKey(string BasisName, int AtomNumber) : IEquatable<BasisKey>
    {
        public readonly bool Equals(BasisKey other)
        {
            return BasisName == other.BasisName && AtomNumber == other.AtomNumber;
        }

        public override readonly int GetHashCode()
        {
            return HashCode.Combine(AtomNumber, BasisName);
        }
    }
}
