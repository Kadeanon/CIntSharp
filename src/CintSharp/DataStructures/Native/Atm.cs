using System;
using System.Collections.Generic;
using System.Diagnostics.Metrics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DataStructures.Native
{
    /// <summary>
    /// Atoms的结构体，存储必要的原子信息以提供给libcint库使用
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct Atm(int charge)
    {

        /// <summary>
        /// Atom charge
        /// </summary>
        public int charge = charge;

        /// <summary>
        /// Ptr of nuclear coordinates, 
        /// storing the offset of the nuclear coordinates in env
        /// </summary>
        public int pointerForCoords = 0;

        /// <summary>
        /// Nuclear model
        /// </summary>
        public NuclearModelKind NuclearModel = NuclearModelKind.PointNuclear;

        /// <summary>
        /// Ptr of ζ, 
        /// is used only when <see cref="NuclearModel"/> 
        /// is <see cref="NuclearModelKind.GaussianNuclear"/> atoms, 
        /// storing the offset of ζ in env
        /// </summary>
        public int pointerForZeta = 0;

        /// <summary>
        /// Unused
        /// </summary>
        public int pointerForFracCharge = 0;

        /// <summary>
        /// Unused
        /// </summary>
        public int reserved = 0;

        public void Reset(int atomNumber = 0)
        {
            charge = atomNumber;
            pointerForCoords = 0;
            NuclearModel = NuclearModelKind.PointNuclear;
            pointerForZeta = 0;
            pointerForFracCharge = 0;
            reserved = 0;
        }

        public readonly int Size => Marshal.SizeOf(this);
    }

    public enum NuclearModelKind : int
    {
        PointNuclear = 1,
        GaussianNuclear = 2
    }
}
