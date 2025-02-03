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
    public class BasisDict(IBasisParser parser, List<double> rawList)
    {

        internal Dictionary<BasisKey, Shells> rawDict = [];

        internal IBasisParser dataParser = parser;

        internal List<double> rawList = rawList;

        public Shells this[BasisKey key]
        {
            get => Get(key);
            set => Set(key, value);
        }

        public bool TryAdd(BasisKey key, [NotNullWhen(true)] out Shells? entry) 
        {
            if (rawDict.TryGetValue(key, out entry))
            {
                return true;
            }
            try
            {
                entry = AddUnChecked(key);
                return true;
            }
            catch (ArgumentException) 
            {
                entry = null;
                return false; 
            }
        }

        public Shells Add(BasisKey key)
        {
            if (rawDict.TryGetValue(key, out Shells? entry))
            {
                return entry;
            }
            try
            {
                return AddUnChecked(key);
            }
            catch (ArgumentException)
            {
                throw new ArgumentException
                    ($"Failed to add basis {key.BasisName} for atom {key.AtomNumber}");
            }
        }

        public bool ContainsKey(BasisKey key)
        {
            return rawDict.ContainsKey(key);
        }

        public bool ContainsKey(string basisName, int atomNumber)
        {
            return rawDict.ContainsKey(new(basisName, atomNumber));
        }

        Shells Get(BasisKey key)
        {
            if (rawDict.TryGetValue(key, out Shells? entry))
            {
                return entry;
            }
            try 
            {
                entry = dataParser.ParseBasis(key, rawList);
                rawDict.Add(key, entry);
                return entry;
            }
            catch(Exception e) 
            {
                throw new ArgumentException($"Failed to get basis {key.BasisName} for atom {key.AtomNumber}", e);
            }
        }

        void Set(BasisKey key, Shells entry)
        {
            rawDict[key] = entry;
        }

        Shells AddUnChecked(BasisKey key)
        {
            var entry = dataParser.ParseBasis(key, rawList);
            rawDict.Add(key, entry);
            return entry;
        }


        public readonly struct BasisKey(string basisName, int atomNumber) : IEquatable<BasisKey>
        {
            public string BasisName { get; } = basisName;

            public int AtomNumber { get; } = atomNumber;

            public readonly bool Equals(BasisKey other)
            {
                return BasisName == other.BasisName && AtomNumber == other.AtomNumber;
            }

            public override readonly bool Equals(object? obj)
            {
                return obj is BasisKey other && Equals(other);
            }

            public override readonly int GetHashCode()
            {
                return HashCode.Combine(AtomNumber, BasisName);
            }

            public static bool operator ==(BasisKey left, BasisKey right)
            {
                return left.Equals(right);
            }

            public static bool operator !=(BasisKey left, BasisKey right)
            {
                return !left.Equals(right);
            }

            public static implicit operator BasisKey(Atom atom)
            {
                return new(atom.BasisName, atom.AtomNumber);
            }
        }
    }
}
