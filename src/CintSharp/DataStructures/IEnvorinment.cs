using CintSharp.DataStructures.Native;
using System.Runtime.CompilerServices;

namespace CintSharp.DataStructures
{
    public interface IEnvorinment
    {
        public Atm[] Atms { get; }

        public Bas[] Bases { get; }

        public double[] Envs { get; }

        public ref EnvHeader EnvHeader => ref EnvHeader.FromSpan(Envs);

        public int Natm { get; }

        public int Nbas { get; }

        public Vector3 AtomCoord(int atomIndex)
        {
            return Unsafe.As<double, Vector3>(ref Envs[Atms[atomIndex].pointerForCoords]);
        }

        public HeaderRinvScope RinvAt(int atomIndex)
        {
            return new(this, atomIndex);
        }

        public HeaderRinvScope RinvAt(Vector3 rinvOrigin)
        {
            return new(this, rinvOrigin);
        }
    }
}
