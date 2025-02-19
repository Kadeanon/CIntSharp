using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace CintSharp.DataStructures
{
    public struct Vector3(double x, double y, double z):
        IEquatable<Vector3>,
        IAdditionOperators<Vector3, Vector3, Vector3>,
        IAdditiveIdentity<Vector3, Vector3>,
        ISubtractionOperators<Vector3, Vector3, Vector3>,
        IUnaryNegationOperators<Vector3, Vector3>,
        IMultiplyOperators<Vector3, Vector3, double>,
        IMultiplyOperators<Vector3, double, Vector3>,
        IDivisionOperators<Vector3, double, Vector3>
    {
        public double X { get; set; } = x;

        public double Y { get; set; } = y;

        public double Z { get; set; } = z;

        public static Vector3 Zero { get; } = new(0, 0, 0);

        public double Length => Math.Sqrt(this * this);

        public double LengthSq => this * this;

        public Span<double> AsSpan() => MemoryMarshal.CreateSpan(ref Unsafe.As<Vector3, double>(ref this), 3);

        public static Vector3 AdditiveIdentity => Zero;

        public static Vector3 operator +(Vector3 a, Vector3 b)
        {
            return new Vector3(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        }

        public static Vector3 operator -(Vector3 a, Vector3 b)
        {
            return new Vector3(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }

        public static Vector3 operator -(Vector3 a)
        {
            return new Vector3(-a.X, -a.Y, -a.Z);
        }

        public static Vector3 operator *(Vector3 a, double b)
        {
            return new Vector3(a.X * b, a.Y * b, a.Z * b);
        }

        public static double operator *(Vector3 a, Vector3 b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        }

        public static Vector3 operator /(Vector3 a, double b)
        {
            return new Vector3(a.X / b, a.Y / b, a.Z / b);
        }

        public readonly bool Equals(Vector3 other)
        {
            return X == other.X && Y == other.Y && Z == other.Z;
        }
    }
}
