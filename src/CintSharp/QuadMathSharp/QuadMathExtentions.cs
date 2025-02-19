
using System.Runtime.InteropServices;
using __complex128 = QuadMathSharp.Complex128;
using System.Text;

namespace QuadMathSharp;

internal static class QuadMathExtentions
{
    /*
     * 
extern "C" {

	// Prototypes for real functions
    extern __float128 acosq(__float128);
	extern __float128 acoshq(__float128);
	extern __float128 asinq(__float128);
	extern __float128 asinhq(__float128);
	extern __float128 atanq(__float128);
	extern __float128 atanhq(__float128);
	extern __float128 atan2q(__float128, __float128);
	extern __float128 cbrtq(__float128);
	extern __float128 ceilq(__float128);
	extern __float128 copysignq(__float128, __float128);
	extern __float128 coshq(__float128);
	extern __float128 cosq(__float128);
	extern __float128 erfq(__float128);
	extern __float128 erfcq(__float128);
	extern __float128 exp2q(__float128);
	extern __float128 expq(__float128);
	extern __float128 expm1q(__float128);
	extern __float128 fabsq(__float128);
	extern __float128 fdimq(__float128, __float128);
	extern int finiteq(__float128);
	extern __float128 floorq(__float128);
	extern __float128 fmaq(__float128, __float128, __float128);
	extern __float128 fmaxq(__float128, __float128);
	extern __float128 fminq(__float128, __float128);
	extern __float128 fmodq(__float128, __float128);
	extern __float128 frexpq(__float128, int*);
	extern __float128 hypotq(__float128, __float128);
	extern int isinfq(__float128);
	extern int ilogbq(__float128);
	extern int isnanq(__float128);
	extern int issignalingq(__float128);
	extern __float128 j0q(__float128);
	extern __float128 j1q(__float128);
	extern __float128 jnq(int, __float128);
	extern __float128 ldexpq(__float128, int);
	extern __float128 lgammaq(__float128);
	extern long long int llrintq(__float128);
	extern long long int llroundq(__float128);
	extern __float128 logbq(__float128);
	extern __float128 logq(__float128);
	extern __float128 log10q(__float128);
	extern __float128 log2q(__float128);
	extern __float128 log1pq(__float128);
	extern long int lrintq(__float128);
	extern long int lroundq(__float128);
	extern __float128 modfq(__float128, __float128*);
	extern __float128 nanq(const char*);
	extern __float128 nearbyintq(__float128);
	extern __float128 nextafterq(__float128, __float128);
	extern __float128 powq(__float128, __float128);
	extern __float128 remainderq(__float128, __float128);
	extern __float128 remquoq(__float128, __float128, int*);
	extern __float128 rintq(__float128);
	extern __float128 roundq(__float128);
	extern __float128 scalblnq(__float128, long int);
	extern __float128 scalbnq(__float128, int);
	extern int signbitq(__float128);
	extern void sincosq(__float128, __float128*, __float128*);
	extern __float128 sinhq(__float128);
	extern __float128 sinq(__float128);
	extern __float128 sqrtq(__float128);
	extern __float128 tanq(__float128);
	extern __float128 tanhq(__float128);
	extern __float128 tgammaq(__float128);
	extern __float128 truncq(__float128);
	extern __float128 y0q(__float128);
	extern __float128 y1q(__float128);
	extern __float128 ynq(int, __float128);


	// Prototypes for complex functions
	extern __float128 cabsq(__complex128);
	extern __float128 cargq(__complex128);
	extern __float128 cimagq(__complex128);
	extern __float128 crealq(__complex128);
	extern __complex128 cacosq(__complex128);
	extern __complex128 cacoshq(__complex128);
	extern __complex128 casinq(__complex128);
	extern __complex128 casinhq(__complex128);
	extern __complex128 catanq(__complex128);
	extern __complex128 catanhq(__complex128);
	extern __complex128 ccosq(__complex128);
	extern __complex128 ccoshq(__complex128);
	extern __complex128 cexpq(__complex128);
	extern __complex128 cexpiq(__float128);
	extern __complex128 clogq(__complex128);
	extern __complex128 clog10q(__complex128);
	extern __complex128 conjq(__complex128);
	extern __complex128 cpowq(__complex128, __complex128);
	extern __complex128 cprojq(__complex128);
	extern __complex128 csinq(__complex128);
	extern __complex128 csinhq(__complex128);
	extern __complex128 csqrtq(__complex128);
	extern __complex128 ctanq(__complex128);
	extern __complex128 ctanhq(__complex128);


	// Prototypes for string <-> __float128 conversion functions
	extern __float128 strtoflt128(const char*, char**);
	extern int quadmath_snprintf(char* str, size_t size,

        const char* format, ...);


	// Macros 
#define FLT128_MAX 1.18973149535723176508575932662800702e4932Q
#define FLT128_MIN 3.36210314311209350626267781732175260e-4932Q
#define FLT128_EPSILON 1.92592994438723585305597794258492732e-34Q
#define FLT128_DENORM_MIN 6.475175119438025110924438958227646552e-4966Q
#define FLT128_MANT_DIG 113
#define FLT128_MIN_EXP (-16381)
#define FLT128_MAX_EXP 16384
#define FLT128_DIG 33
#define FLT128_MIN_10_EXP (-4931)
#define FLT128_MAX_10_EXP 4932


#define HUGE_VALQ __builtin_huge_valq()
// The following alternative is valid, but brings the warning:
// (floating constant exceeds range of ¡®__float128¡¯)  
// #define HUGE_VALQ (__extension__ 0x1.0p32767Q) 

#define M_Eq		2.718281828459045235360287471352662498Q   e 
#define M_LOG2Eq	1.442695040888963407359924681001892137Q   log_2 e 
#define M_LOG10Eq	0.434294481903251827651128918916605082Q   log_10 e 
#define M_LN2q		0.693147180559945309417232121458176568Q   log_e 2 
#define M_LN10q		2.302585092994045684017991454684364208Q   log_e 10 
#define M_PIq		3.141592653589793238462643383279502884Q   pi 
#define M_PI_2q		1.570796326794896619231321691639751442Q   pi/2 
#define M_PI_4q		0.785398163397448309615660845819875721Q   pi/4 
#define M_1_PIq		0.318309886183790671537767526745028724Q   1/pi 
#define M_2_PIq		0.636619772367581343075535053490057448Q   2/pi 
#define M_2_SQRTPIq	1.128379167095512573896158903121545172Q   2/sqrt(pi) 
#define M_SQRT2q	1.414213562373095048801688724209698079Q   sqrt(2) 
#define M_SQRT1_2q	0.707106781186547524400844362104849039Q   1/sqrt(2) 
}
	*/

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 acosq(__float128 x);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 acoshq(__float128 x);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 asinq(__float128 x);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 asinhq(__float128 x);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 atanq(__float128 x);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 atanhq(__float128 x);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 atan2q(__float128 x, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 cbrtq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 ceilq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 copysignq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 coshq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 cosq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 erfq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 erfcq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 exp2q(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 expq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 expm1q(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 fabsq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 fdimq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int finiteq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 floorq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 fmaq(__float128 x, __float128 y, __float128 z);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 fmaxq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 fminq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 fmodq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 frexpq(__float128 x, ref int exp);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 hypotq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int isinfq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int ilogbq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int isnanq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int issignalingq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 j0q(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 j1q(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 jnq(int n, __float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 ldexpq(__float128 x, int exponent);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 lgammaq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern long llrintq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern long llroundq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 logbq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 logq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 log10q(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 log2q(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 log1pq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int lrintq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int lroundq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 modfq(__float128 x, ref __float128 integer);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 nanq(ref readonly char str);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 nearbyintq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 nextafterq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 powq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 remainderq(__float128 left, __float128 right);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 remquoq(__float128 x, __float128 y, ref int quo);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 rintq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 roundq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 scalblnq(__float128 arg, int exp);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 scalbnq(__float128 arg, int exp);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int signbitq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern void sincosq(__float128 x, ref __float128 sin, ref __float128 cos);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 sinhq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 sinq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 sqrtq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 tanq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 tanhq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 tgammaq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 truncq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 y0q(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 y1q(__float128 x);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 ynq(int n, __float128 x);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 cabsq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 cargq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 cimagq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 crealq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 cacosq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 cacoshq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 casinq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 casinhq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 catanq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 catanhq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 ccosq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 ccoshq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 cexpq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 cexpiq(__float128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 clogq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 clog10q(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 conjq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 cpowq(__complex128 x, __complex128 m);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 cprojq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 csinq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 csinhq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 csqrtq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 ctanq(__complex128 x);
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __complex128 ctanhq(__complex128 x);

    #region 
    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl)]
    internal static extern __float128 strtoflt128([MarshalAs(UnmanagedType.LPStr)] string str, IntPtr endPtr);

    [DllImport("Native/libquadmath-0.dll", CallingConvention = CallingConvention.Cdecl, SetLastError = true)]
    internal static extern int quadmath_snprintf([MarshalAs(UnmanagedType.LPStr)] StringBuilder str, ulong size,
        [MarshalAs(UnmanagedType.LPStr)] string format, [In] __float128  arg);


    public const int FLT128_MANT_DIG = 113;
    public const int FLT128_MIN_EXP = (-16381);
    public const int FLT128_MAX_EXP = 16384;
    public const int FLT128_DIG = 33;
    public const int FLT128_MIN_10_EXP = -4931;
    public const int FLT128_MAX_10_EXP = 4932;

    public static readonly __float128 FLT128_MAX = __float128.FromString("1.18973149535723176508575932662800702e4932");
    public static readonly __float128 FLT128_MIN = __float128.FromString("3.36210314311209350626267781732175260e-4932");
    public static readonly __float128 FLT128_EPSILON = __float128.FromString("1.92592994438723585305597794258492732e-34");
    public static readonly __float128 FLT128_DENORM_MIN = __float128.FromString("6.475175119438025110924438958227646552e-4966");

    /// <summary>
    /// Euler's number
    /// </summary>
    public static readonly __float128 M_Eq = __float128.FromString("2.718281828459045235360287471352662498");

    /// <summary>
    /// log_2 e
    /// </summary>
    public static readonly __float128 M_LOG2Eq = __float128.FromString("1.442695040888963407359924681001892137");

    /// <summary>
    /// log_10 e
    /// </summary>
    public static readonly __float128 M_LOG10Eq = __float128.FromString("0.434294481903251827651128918916605082");

    /// <summary>
    /// log_e 2
    /// </summary>
    public static readonly __float128 M_LN2q = __float128.FromString("0.693147180559945309417232121458176568");

    /// <summary>
    /// log_e 10
    /// </summary>
    public static readonly __float128 M_LN10q = __float128.FromString("2.302585092994045684017991454684364208");

    /// <summary>
    /// pi
    /// </summary>
    public static readonly __float128 M_PIq = __float128.FromString("3.141592653589793238462643383279502884");

    /// <summary>
    /// pi/2
    /// </summary>
    public static readonly __float128 M_PI_2q = __float128.FromString("1.570796326794896619231321691639751442");

    /// <summary>
    /// pi/4
    /// </summary>
    public static readonly __float128 M_PI_4q = __float128.FromString("0.785398163397448309615660845819875721");

    /// <summary>
    /// 1/pi
    /// </summary>
    public static readonly __float128 M_1_PIq = __float128.FromString("0.318309886183790671537767526745028724");

    /// <summary>
    /// 2/pi
    /// </summary>
    public static readonly __float128 M_2_PIq = __float128.FromString("0.636619772367581343075535053490057448");

    /// <summary>
    /// 2/sqrt(pi)
    /// </summary>
    public static readonly __float128 M_2_SQRTPIq = __float128.FromString("1.128379167095512573896158903121545172");

    /// <summary>
    /// sqrt(2)
    /// </summary>
    public static readonly __float128 M_SQRT2q = __float128.FromString("1.414213562373095048801688724209698079");

    /// <summary>
    /// 1/sqrt(2)
    /// </summary>
    public static readonly __float128 M_SQRT1_2q = __float128.FromString("0.707106781186547524400844362104849039");
    #endregion
}