using System;
using System.Buffers;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Linq;
using System.Runtime.InteropServices.Marshalling;
using System.Text;
using System.Threading.Tasks;

namespace QuadMathSharp
{
    public partial struct Float128
    {
        public static __float128 Parse(string s) 
            => Parse(s.AsSpan(), NumberStyles.Float | NumberStyles.AllowThousands, provider: null);

        public static __float128 Parse(string s, NumberStyles style)
            => Parse(s.AsSpan(), style, provider: null);

        public static __float128 Parse(string s, IFormatProvider? provider) 
            => Parse(s.AsSpan(), NumberStyles.Float | NumberStyles.AllowThousands, provider);

        public static __float128 Parse(string s, NumberStyles style, IFormatProvider? provider) 
            => Parse(s.AsSpan(), style, provider);

        public static __float128 Parse(ReadOnlySpan<char> s, IFormatProvider? provider) 
            => Parse(s, NumberStyles.Float | NumberStyles.AllowThousands, provider);

        public static __float128 Parse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider)
        {
            if (!TryParse(s, style, provider, out var result))
            {
                throw new FormatException($"Bad format: cannot parse Float128 struct with {s.ToString()} format.");
            }
            return result;
        }

        public override readonly string ToString()
        {
            return ToString(null, NumberFormatInfo.CurrentInfo);
        }

        public readonly string ToString(string? format, IFormatProvider? formatProvider)
        {
            var info = NumberFormatInfo.GetInstance(formatProvider);
            const int bufferSize = 1024;
            char c = ParseFormatSpecifier(format, out int digits);
            if (c == '\0') c = format[0];
            c = Char.ToLower(c);
            string fmt;
            switch (c)
            {
                case 'c':
                    throw new NotSupportedException();
                case 'e':
                    if (digits == -1) digits = 6;
                    fmt = $"%.{digits}Qe";
                    break;
                case 'f':
                    if (digits == -1) digits = info.NumberDecimalDigits;
                    fmt = $"%.{digits}Qf";
                    break;
                case 'g':
                    if (digits == -1) digits = QuadMathExtentions.FLT128_DIG;
                    fmt = $"%.{digits}Qe";
                    break;
                case 'r':
                case 'n':
                case 'p':
                    throw new NotImplementedException();
                default:
                    throw new FormatException($"Bad format: '{c}' is not implemented for Float128!");
            }
            using(StringBuilderPool.ScopedBorrow(bufferSize, out var sb))
            {
                int cosumed = QuadMathExtentions.quadmath_snprintf(sb, bufferSize, fmt, this);
                if (cosumed == -1)
                {
                    throw new FormatException($"Bad format: cannot format Float128 struct with {format} format.");
                }
                return sb.ToString();
            }
        }

        internal static char ParseFormatSpecifier(ReadOnlySpan<char> format, out int digits)
        {
            char c = default;
            if (format.Length > 0)
            {
                // If the format begins with a symbol, see if it's a standard format
                // with or without a specified number of digits.
                c = format[0];
                if (char.IsAsciiLetter(c))
                {
                    // Fast path for sole symbol, e.g. "D"
                    if (format.Length == 1)
                    {
                        digits = -1;
                        return c;
                    }

                    if (format.Length == 2)
                    {
                        // Fast path for symbol and single digit, e.g. "X4"
                        int d = format[1] - '0';
                        if ((uint)d < 10)
                        {
                            digits = d;
                            return c;
                        }
                    }
                    else if (format.Length == 3)
                    {
                        // Fast path for symbol and double digit, e.g. "F12"
                        int d1 = format[1] - '0', d2 = format[2] - '0';
                        if ((uint)d1 < 10 && (uint)d2 < 10)
                        {
                            digits = d1 * 10 + d2;
                            return c;
                        }
                    }

                    // Fallback for symbol and any length digits.  The digits value must be >= 0 && <= 999_999_999,
                    // but it can begin with any number of 0s, and thus we may need to check more than 9
                    // digits.  Further, for compat, we need to stop when we hit a null char.
                    int n = 0;
                    int i = 1;
                    while ((uint)i < (uint)format.Length && char.IsAsciiDigit(format[i]))
                    {
                        // Check if we are about to overflow past our limit of 9 digits
                        if (n >= 100_000_000)
                        {
                            throw new FormatException($"Bad format ");
                        }
                        n = (n * 10) + format[i++] - '0';
                    }

                    // If we're at the end of the digits rather than having stopped because we hit something
                    // other than a digit or overflowed, return the standard format info.
                    if ((uint)i >= (uint)format.Length || format[i] == '\0')
                    {
                        digits = n;
                        return c;
                    }
                }
            }

            // Default empty format to be "G"; custom format is signified with '\0'.
            digits = -1;
            return format.Length == 0 || c == '\0' ? // For compat, treat '\0' as the end of the specifier, even if the specifier extends beyond it.
                'G' :
                '\0';
        }

        public readonly bool TryFormat(Span<char> destination, out int charsWritten, ReadOnlySpan<char> format, IFormatProvider? provider)
        {
            var info = NumberFormatInfo.GetInstance(provider);
            int bufferSize = destination.Length;
            char c = ParseFormatSpecifier(format, out int digits);
            if (c == '\0') c = format[0];
            c = Char.ToLower(c);
            string fmt;
            switch (c)
            {
                case 'c':
                    throw new NotSupportedException();
                case 'e':
                    if (digits == -1) digits = 6;
                    fmt = $"%.{digits}Qe";
                    break;
                case 'f':
                    if (digits == -1) digits = info.NumberDecimalDigits;
                    fmt = $"%.{digits}Qf";
                    break;
                case 'g':
                    if (digits == -1) digits = QuadMathExtentions.FLT128_DIG;
                    fmt = $"%.{digits}Qe";
                    break;
                case 'r':
                case 'n':
                case 'p':
                    throw new NotImplementedException();
                default:
                    throw new FormatException($"Bad format: '{c}' is not implemented for Float128!");
            }
            using (StringBuilderPool.ScopedBorrow(bufferSize, out var sb))
            {
                charsWritten = QuadMathExtentions.quadmath_snprintf(sb, (ulong)bufferSize, fmt, this);
                if (charsWritten == -1)
                {
                    charsWritten = 0;
                    return false;
                }
                else 
                {
                    var span = sb.ToString().AsSpan();
                    span[..charsWritten].CopyTo(destination);
                    return true;
                }
            }
        }

        public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, [MaybeNullWhen(false)] out __float128 result)
        {
            throw new NotImplementedException();
        }

        public static bool TryParse([NotNullWhen(true)] string? s, NumberStyles style, IFormatProvider? provider, [MaybeNullWhen(false)] out __float128 result)
            => TryParse(s.AsSpan(), style, provider, out result);

        public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, [MaybeNullWhen(false)] out __float128 result)
            => TryParse(s, NumberStyles.Float | NumberStyles.AllowThousands, provider, out result);

        public static bool TryParse([NotNullWhen(true)] string? s, IFormatProvider? provider, [MaybeNullWhen(false)] out __float128 result)
            => TryParse(s.AsSpan(), NumberStyles.Float | NumberStyles.AllowThousands, provider, out result);

    }
}
