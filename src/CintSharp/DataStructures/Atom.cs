using CintSharp.DataStructures.Native;
using NCDK;
using NCDK.Default;
using NCDK.Numerics;
using NCDK.Tools;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DataStructures
{
    public struct Atom
    {

        public Vector3 position;

        public ChemicalElement Element { get; }

        public string BasisName { get; set; }

        public double X 
        { 
            readonly get => position.X; 
            set => position.X = value; 
        }

        public double Y
        {
            readonly get => position.Y;
            set => position.Y = value;
        }

        public double Z
        {
            readonly get => position.Z;
            set => position.Z = value;
        }

        public readonly int AtomNumber => Element.AtomicNumber;

        public Atom(ChemicalElement element, double x, double y, double z, string basisName = "")
        {
            Element = element;
            X = x;
            Y = y;
            Z = z;
            BasisName = basisName;
        }

        public Atom(int atomNumber, double x, double y, double z, string basisName = "")
            : this(AtomicNumberToElement(atomNumber), x, y, z, basisName) { }

        public Atom(string atomSymbol, double x, double y, double z, string basisName = "")
            : this(SymbolToElement(atomSymbol), x, y, z, basisName) { }

        public Atom(IAtom atom, string basisName = "")
        {
            Element = SymbolToElement(atom.Symbol);
            Vector3? atomPosition = atom.Point3D;
            if(!atomPosition.HasValue)
            {
                throw new ArgumentException("Atom does not have a 3d position");
            }
            position = atomPosition.Value;
            BasisName = basisName;
        }

        public void CopyToAtm(ref Atm atm, List<double> envs)
        {
            atm.Reset(AtomNumber);
            atm.pointerForCoords = envs.Count;
            envs.AddRange(MemoryMarshal.CreateSpan(ref Unsafe.As<Vector3, double>(ref position), 3));
        }

        public static ChemicalElement SymbolToElement(string symbol)
        {
            return ChemicalElement.OfSymbol(symbol);
        }

        public static ChemicalElement AtomicNumberToElement(int atomicNumber)
        {
            return ChemicalElement.Of(atomicNumber);
        }
    }
}
