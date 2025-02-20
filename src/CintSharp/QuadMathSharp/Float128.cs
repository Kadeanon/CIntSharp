global using __float128 = QuadMathSharp.Float128;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;


namespace QuadMathSharp
{
    [StructLayout(LayoutKind.Explicit, Size = 16)]
    public partial struct Float128 : IBinaryFloatingPointIeee754<__float128>
    {
        [FieldOffset(0)]
        internal double low;

        [FieldOffset(8)]
        internal double high;

        internal ref Int128 int128Value => ref ByteMarshal<Int128>(0);

        internal ref uint mantissa3 => ref ByteMarshal<uint>(0);

        internal ref uint mantissa2 => ref ByteMarshal<uint>(4);

        internal ref uint mantissa1 => ref ByteMarshal<uint>(8);

        internal ref ushort mantissa0 => ref ByteMarshal<ushort>(12);

        internal ref short exponent => ref ByteMarshal<short>(14);

        internal bool Sign => exponent < 0;

        internal bool QuietNaN => (mantissa0 & 0x8000) != 0;

        private Float128(Int128 value)
        {
            int128Value = value;
        }

        public static implicit operator Float128(double value) => SoftFpExtentions.__extenddftf2(value);

        public static implicit operator Float128(float value) => SoftFpExtentions.__extendsftf2(value);

        public static explicit operator double(Float128 value) => SoftFpExtentions.__trunctfdf2(value);

        public static explicit operator float(Float128 value) => SoftFpExtentions.__trunctfsf2(value);

        internal ref T ByteMarshal<T>(int index) where T : unmanaged
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(index, 0, nameof(index));
            ArgumentOutOfRangeException.ThrowIfGreaterThan(index, 15, nameof(index));
            Span<byte> span = MemoryMarshal.CreateSpan(ref Unsafe.As<__float128, byte>(ref this), 16);
            return ref Unsafe.As<byte, T>(ref span[index]);
        }

        public static __float128 Two { get; } = FromString("2");

        public static __float128 Ten { get; } = FromString("10");

        public static __float128 Epsilon => QuadMathExtentions.FLT128_EPSILON;

        public static __float128 NaN => throw new NotImplementedException();

        public static __float128 NegativeInfinity => throw new NotImplementedException();

        public static __float128 NegativeZero => throw new NotImplementedException();

        public static __float128 PositiveInfinity => throw new NotImplementedException();

        public static __float128 NegativeOne => throw new NotImplementedException();

        public static __float128 E => QuadMathExtentions.M_Eq;

        public static __float128 Pi => QuadMathExtentions.M_PIq;

        public static __float128 Tau => throw new NotImplementedException();

        public static __float128 One { get; } = FromString("1");

        public static int Radix => 2;

        public static __float128 Zero { get; } = FromString("0");

        public static __float128 AdditiveIdentity => Zero;

        public static __float128 MultiplicativeIdentity => One;

        public static __float128 Abs(__float128 value)
        {
            return QuadMathExtentions.fabsq(value);
        }

        public static __float128 Acos(__float128 x)
        {
            return QuadMathExtentions.acosq(x);
        }

        public static __float128 Acosh(__float128 x)
        {
            return QuadMathExtentions.acoshq(x);
        }

        public static __float128 AcosPi(__float128 x)
        {
            return QuadMathExtentions.acosq(x) / Pi;
        }

        public static __float128 Asin(__float128 x)
        {
            return QuadMathExtentions.asinq(x);
        }

        public static __float128 Asinh(__float128 x)
        {
            return QuadMathExtentions.asinhq(x);
        }

        public static __float128 AsinPi(__float128 x)
        {
            return Asin(x) / Pi;
        }

        public static __float128 Atan(__float128 x)
        {
            return QuadMathExtentions.atanq(x);
        }

        public static __float128 Atan2(__float128 y, __float128 x)
        {
            return QuadMathExtentions.atan2q(y, x);
        }

        public static __float128 Atan2Pi(__float128 y, __float128 x)
        {
            return Atan2(y, x) / Pi;
        }

        public static __float128 Atanh(__float128 x)
        {
            return QuadMathExtentions.atanhq(x);
        }

        public static __float128 AtanPi(__float128 x)
        {
            return Atan(x) / Pi;
        }

        public static __float128 BitDecrement(__float128 x)
        {
            throw new NotImplementedException();
        }

        public static __float128 BitIncrement(__float128 x)
        {
            throw new NotImplementedException();
        }

        public static __float128 Cbrt(__float128 x)
        {
            return QuadMathExtentions.cbrtq(x);
        }

        public static __float128 Cos(__float128 x)
        {
            return QuadMathExtentions.cosq(x);
        }

        public static __float128 Cosh(__float128 x)
        {
            return QuadMathExtentions.coshq(x);
        }

        public static __float128 CosPi(__float128 x)
        {
            return Cos(x * Pi);
        }

        public static __float128 Erf(__float128 x)
        {
            return QuadMathExtentions.erfq(x);
        }

        public static __float128 Erfc(__float128 x)
        {
            return QuadMathExtentions.erfcq(x);
        }

        public static __float128 Exp(__float128 x)
        {
            return QuadMathExtentions.expq(x);
        }

        public static __float128 Exp10(__float128 x)
        {
            return Exp(x * Log(FromString("10")));
        }

        public static __float128 Exp2(__float128 x)
        {
            return QuadMathExtentions.exp2q(x);
        }

        public static Float128 FromString(string str)
        {
            Float128 value = QuadMathExtentions.strtoflt128(str, (nint)0);
            return value;
        }

        public static __float128 FusedMultiplyAdd(__float128 left, __float128 right, __float128 addend)
        {
            return QuadMathExtentions.fmaq(left, right, addend);
        }

        public static __float128 Hypot(__float128 x, __float128 y)
        {
            throw new NotImplementedException();
        }

        public static __float128 Ieee754Remainder(__float128 left, __float128 right)
        {
            return QuadMathExtentions.remainderq(left, right);
        }

        public static int ILogB(__float128 x)
        {
            return QuadMathExtentions.ilogbq(x);
        }

        public static bool IsCanonical(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsComplexNumber(__float128 value) => false;

        public static bool IsEvenInteger(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsFinite(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsImaginaryNumber(__float128 value) => false;

        public static bool IsInfinity(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsInteger(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsNaN(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsNegative(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsNegativeInfinity(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsNormal(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsOddInteger(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsPositive(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsPositiveInfinity(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsPow2(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsRealNumber(__float128 value) => true;

        public static bool IsSubnormal(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static bool IsZero(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static __float128 Log(__float128 x)
        {
            return QuadMathExtentions.logq(x);
        }

        public static __float128 Log(__float128 x, __float128 newBase)
        {
            return QuadMathExtentions.logq(x) / QuadMathExtentions.logq(newBase);
        }

        public static __float128 Log10(__float128 x)
        {
            throw new NotImplementedException();
        }

        public static __float128 Log2(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static __float128 MaxMagnitude(__float128 x, __float128 y)
        {
            throw new NotImplementedException();
        }

        public static __float128 MaxMagnitudeNumber(__float128 x, __float128 y)
        {
            throw new NotImplementedException();
        }

        public static __float128 MinMagnitude(__float128 x, __float128 y)
        {
            throw new NotImplementedException();
        }

        public static __float128 MinMagnitudeNumber(__float128 x, __float128 y)
        {
            throw new NotImplementedException();
        }

        public static __float128 Pow(__float128 x, __float128 y)
        {
            throw new NotImplementedException();
        }

        public static __float128 RootN(__float128 x, int n)
        {
            throw new NotImplementedException();
        }

        public static __float128 Round(__float128 x, int digits, MidpointRounding mode)
        {
            throw new NotImplementedException();
        }

        public static __float128 ScaleB(__float128 x, int n)
        {
            return QuadMathExtentions.scalbnq(x, n);
        }

        public static __float128 Sin(__float128 x)
        {
            return QuadMathExtentions.sinq(x);
        }

        public static (__float128 Sin, __float128 Cos) SinCos(__float128 x)
        {
            __float128 sin = default;
            __float128 cos = default;
            QuadMathExtentions.sincosq(x, ref sin, ref cos);
            return (sin, cos);
        }

        public static (__float128 SinPi, __float128 CosPi) SinCosPi(__float128 x)
        {
            __float128 sin = default;
            __float128 cos = default;
            QuadMathExtentions.sincosq(x * Pi, ref sin, ref cos);
            return (sin, cos);
        }

        public static __float128 Sinh(__float128 x)
        {
            return QuadMathExtentions.sinhq(x);
        }

        public static __float128 SinPi(__float128 x)
        {
            return Sin(x * Pi);
        }

        public static __float128 Sqrt(__float128 x)
        {
            return QuadMathExtentions.sqrtq(x);
        }

        public static __float128 Tan(__float128 x)
        {
            return QuadMathExtentions.tanq(x);
        }

        public static __float128 Tanh(__float128 x)
        {
            return QuadMathExtentions.tanhq(x);
        }

        public static __float128 TanPi(__float128 x)
        {
            return Tan(x * Pi);
        }

        static bool INumberBase<__float128>.TryConvertFromChecked<TOther>(TOther value, out __float128 result)
        {
            throw new NotImplementedException();
        }

        static bool INumberBase<__float128>.TryConvertFromSaturating<TOther>(TOther value, out __float128 result)
        {
            throw new NotImplementedException();
        }

        static bool INumberBase<__float128>.TryConvertFromTruncating<TOther>(TOther value, out __float128 result)
        {
            throw new NotImplementedException();
        }

        static bool INumberBase<__float128>.TryConvertToChecked<TOther>(__float128 value, out TOther result)
        {
            throw new NotImplementedException();
        }

        static bool INumberBase<__float128>.TryConvertToSaturating<TOther>(__float128 value, out TOther result)
        {
            throw new NotImplementedException();
        }

        static bool INumberBase<__float128>.TryConvertToTruncating<TOther>(__float128 value, out TOther result)
        {
            throw new NotImplementedException();
        }

        public int CompareTo(object? obj)
        {
            if (obj is __float128 other)
            {
                return CompareTo(other);
            }
            else
            {
                ArgumentException.ThrowIfNullOrEmpty(null, nameof(obj));
                return -1;
            }
        }

        public int CompareTo(__float128 other)
        {
            throw new NotImplementedException();
        }

        public bool Equals(__float128 other)
        {
            return int128Value == other.int128Value;
        }

        public int GetExponentByteCount()
        {
            throw new NotImplementedException();
        }

        public int GetExponentShortestBitLength()
        {
            throw new NotImplementedException();
        }

        public int GetSignificandBitLength()
        {
            throw new NotImplementedException();
        }

        public int GetSignificandByteCount()
        {
            throw new NotImplementedException();
        }

        public bool TryWriteExponentBigEndian(Span<byte> destination, out int bytesWritten)
        {
            throw new NotImplementedException();
        }

        public bool TryWriteExponentLittleEndian(Span<byte> destination, out int bytesWritten)
        {
            throw new NotImplementedException();
        }

        public bool TryWriteSignificandBigEndian(Span<byte> destination, out int bytesWritten)
        {
            throw new NotImplementedException();
        }

        public bool TryWriteSignificandLittleEndian(Span<byte> destination, out int bytesWritten)
        {
            throw new NotImplementedException();
        }

        public static __float128 operator +(__float128 value)
        {
            return value;
        }

        public static __float128 operator +(__float128 left, __float128 right)
        {
            return SoftFpExtentions.__addtf3(left, right);
        }

        public static __float128 operator -(__float128 value)
        {
            return SoftFpExtentions.__negtf2(value);
        }

        public static __float128 operator -(__float128 left, __float128 right)
        {
            return SoftFpExtentions.__subtf3(left, right);
        }

        public static __float128 operator ~(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static __float128 operator ++(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static __float128 operator --(__float128 value)
        {
            throw new NotImplementedException();
        }

        public static __float128 operator *(__float128 left, __float128 right)
        {
            return SoftFpExtentions.__multf3(left, right);
        }

        public static __float128 operator /(__float128 left, __float128 right)
        {
            return SoftFpExtentions.__divtf3(left, right);
        }

        public static __float128 operator %(__float128 left, __float128 right)
        {
            throw new NotImplementedException();
        }

        public static __float128 operator &(__float128 left, __float128 right)
        {
            throw new NotImplementedException();
        }

        public static __float128 operator |(__float128 left, __float128 right)
        {
            throw new NotImplementedException();
        }

        public static __float128 operator ^(__float128 left, __float128 right)
        {
            throw new NotImplementedException();
        }

        public static bool operator ==(__float128 left, __float128 right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(__float128 left, __float128 right)
        {
            return !left.Equals(right);
        }

        public static bool operator <(__float128 left, __float128 right)
        {
            return SoftFpExtentions.__lttf2(left, right) != 0;
        }

        public static bool operator >(__float128 left, __float128 right)
        {
            return SoftFpExtentions.__gttf2(left, right) != 0;
        }

        public static bool operator <=(__float128 left, __float128 right)
        {
            return SoftFpExtentions.__letf2(left, right) != 0;
        }

        public static bool operator >=(__float128 left, __float128 right)
        {
            return SoftFpExtentions.__getf2(left, right) != 0;
        }

        public override bool Equals(object obj)
        {
            return obj is __float128 other && Equals(other);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(int128Value);
        }
    }
}
