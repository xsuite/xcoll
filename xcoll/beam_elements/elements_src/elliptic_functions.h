#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

#include <math.h>

// ---- Mathematical constants (renamed to avoid conflicts) ----
#define MY_PI       3.14159265358979323846
#define MY_PIO2     1.57079632679489661923
//Cephes uses MACHEP 2^(-53) as a tight bound on relative rounding error in many formulas and stopping criteria.
#define MY_MACHEP   1.11022302462515654042E-16
#define MY_MAXNUM   1.7976931348623158E308


//-------------mconf.h--------------------
/* Define if the `long double' type works.  */
#define HAVE_LONG_DOUBLE 1

/* Define as the return type of signal handlers (int or void).  */
#define RETSIGTYPE void

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* Define if floating point words are bigendian.  */
/* #undef FLOAT_WORDS_BIGENDIAN */

/* The number of bytes in a int.  */
#define SIZEOF_INT 4

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H 1

/* Name of package */
#define PACKAGE "cephes"

/* Version number of package */
#define VERSION "2.7"

/* Constant definitions for math error conditions
 */

#define DOMAIN		1	/* argument domain error */
#define SING		2	/* argument singularity */
#define OVERFLOW	3	/* overflow range error */
#define UNDERFLOW	4	/* underflow range error */
#define TLOSS		5	/* total loss of precision */
#define PLOSS		6	/* partial loss of precision */

#define EDOM		33
#define ERANGE		34
/* Complex numeral.  */
typedef struct
	{
	double r;
	double i;
	} cmplx;

#ifdef HAVE_LONG_DOUBLE
/* Long double complex numeral.  */
typedef struct
	{
	long double r;
	long double i;
	} cmplxl;
#endif


/* Type of computer arithmetic */

/* PDP-11, Pro350, VAX:
 */
/* #define DEC 1 */

/* Intel IEEE, low order words come first:
 */
/* #define IBMPC 1 */

/* Motorola IEEE, high order words come first
 * (Sun 680x0 workstation):
 */
/* #define MIEEE 1 */

/* UNKnown arithmetic, invokes coefficients given in
 * normal decimal format.  Beware of range boundary
 * problems (MACHEP, MAXLOG, etc. in const.c) and
 * roundoff problems in pow.c:
 * (Sun SPARCstation)
 */
#define UNK 1

/* If you define UNK, then be sure to set BIGENDIAN properly. */
#ifdef FLOAT_WORDS_BIGENDIAN
#define BIGENDIAN 1
#else
#define BIGENDIAN 0
#endif
/* Define this `volatile' if your compiler thinks
 * that floating point arithmetic obeys the associative
 * and distributive laws.  It will defeat some optimizations
 * (but probably not enough of them).
 *
 * #define VOLATILE volatile
 */
#define VOLATILE

/* For 12-byte long doubles on an i386, pad a 16-bit short 0
 * to the end of real constants initialized by integer arrays.
 *
 * #define XPD 0,
 *
 * Otherwise, the type is 10 bytes long and XPD should be
 * defined blank (e.g., Microsoft C).
 *
 * #define XPD
 */
#define XPD 0,

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

/* Define to ask for infinity support, else undefine. */
/* #define INFINITIES 1 */

/* Define to ask for support of numbers that are Not-a-Number,
   else undefine.  This may automatically define INFINITIES in some files. */
/* #define NANS 1 */

/* Define to distinguish between -0.0 and +0.0.  */
#define MINUSZERO 1

/* Define 1 for ANSI C atan2() function
   See atan.c and clog.c. */
#define ANSIC 1

/* Get ANSI function prototypes, if you want them. */

// Forward declaration matches our implementation
GPUFUN int mtherr(char *name, int code);
GPUFUN int get_cephes_errno(void);

/* Variable for error reporting.  See mtherr.c.  */
extern int merror;


//--------------mtherr.c----------------

static int cephes_errno = 0;

// Implementation for GPU/CPU
GPUFUN int mtherr(char *name, int code)
{
    (void)name; // unused parameter

    if (code <= 0 || code > PLOSS)
        code = PLOSS;

    cephes_errno = code;
    return 0;
}

GPUFUN int get_cephes_errno(void)
{
    int ret = cephes_errno;
    cephes_errno = 0;
    return ret;
}







// IS THIS POLYNOMIAL EVALUATION ENOUGH? https://lost-contact.mit.edu/afs/lngs.infn.it/experiment/cuore/soft_node101/scipy-0.16.0/scipy/special/cephes/polevl.h
//---------------------------------------
// ---- Polynomial evaluation (Cephes) ----
GPUFUN  double polevl(double x, const double coef[], int N) {
    double ans = coef[0];
    for (int i = 1; i <= N; ++i) {
        ans = ans * x + coef[i];
    }
    return ans;
}

GPUFUN double p1evl(double x, const double coef[], int N) {
    double ans = x + coef[0];
    for (int i = 1; i < N; ++i) {
        ans = ans * x + coef[i];
    }
    return ans;
}

// ELLPK needs fixing????
//-------------------
// ---- ellpk: complete elliptic integral of the first kind ----
// Complete elliptic integral of the first kind K(m) with x = m1 = 1 - m.

// Canonical UNK coefficients (portable IEEE double)
static const double ELLPK_P_UNK[11] = {
    1.37982864606273237150e-4, 2.28025724005875567385e-3,
    7.97404013220415179367e-3, 9.85821379021226008714e-3,
    6.87489687449949877925e-3, 6.18901033637687613229e-3,
    8.79078273952743772254e-3, 1.49380448916805252718e-2,
    3.08851465246711995998e-2, 9.65735902811690126535e-2,
    1.38629436111989062502e0
};

static const double ELLPK_Q_UNK[11] = {
    2.94078955048598507511e-5, 9.14184723865917226571e-4,
    5.94058303753167793257e-3, 1.54850516649762399335e-2,
    2.39089602715924892727e-2, 3.01204715227604046988e-2,
    3.73774314173823228969e-2, 4.88280347570998239232e-2,
    7.03124996963957469739e-2, 1.24999999999870820058e-1,
    4.99999999999999999821e-1
};



#  define ELLPK_P  ELLPK_P_UNK     /* default */
#  define ELLPK_Q  ELLPK_Q_UNK

#define C1 1.3862943611198906188
/* Optional: exact log(4) used in the small-x asymptotic path */

GPUFUN double ellpk(double x)
{
    if (x < 0.0 || x > 1.0) {
        // Domain error, return 0 for safety (Cephes behavior)
        return 0.0;
    }

    if (x > MY_MACHEP) {
        double p = polevl(x, ELLPK_P, 10);
        double q = polevl(x, ELLPK_Q, 10);
        return p - log(x) * q;
    }

    // Very small x: use asymptotic form
    if (x == 0.0) {
        return MY_MAXNUM;  // singularity
    }

    return C1 - 0.5 * log(x);
}


// ---- ellik: incomplete elliptic integral of the first kind ----
// ---- ellik: incomplete elliptic integral of the first kind F(φ | m) ----
GPUFUN double ellik(double phi, double m)
{
    // Trivial case
    if (m == 0.0) {
        return phi;
    }

    double a = 1.0 - m;

    // m → 1: logarithmic singularity
    if (a == 0.0) {
        if (fabs(phi) >= MY_PIO2) {
            return MY_MAXNUM;
        }
        return log(tan((MY_PIO2 + phi) * 0.5));
    }

    // Reduce φ modulo π/2
    int npio2 = (int)floor(phi / MY_PIO2);
    if (npio2 & 1) {
        npio2 += 1;
    }

    double K;
    if (npio2) {
        K   = ellpk(a);
        phi = phi - npio2 * MY_PIO2;
    } else {
        K = 0.0;
    }

    int sign;
    if (phi < 0.0) {
        phi  = -phi;
        sign = -1;
    } else {
        sign = 0;
    }

    double b = sqrt(a);
    double t = tan(phi);

    // Avoid excessive recursion
    if (fabs(t) > 10.0) {
        double e = 1.0 / (b * t);
        if (fabs(e) < 10.0) {
            e = atan(e);
            if (npio2 == 0) {
                K = ellpk(a);
            }
            double temp = K - ellik(e, m);
            double out  = (sign ? -1.0 : 1.0) * temp + npio2 * K;
            return out;
        }
    }

    // AGM iteration
    a = 1.0;
    double c   = sqrt(m);
    int    d   = 1;
    int    mod = 0;

    while (fabs(c / a) > MY_MACHEP) {
        double temp = b / a;
        phi = phi + atan(t * temp) + mod * MY_PI;
        mod = (int)((phi + MY_PIO2) / MY_PI);

        t = t * (1.0 + temp) / (1.0 - temp * t * t);

        c = 0.5 * (a - b);

        double ab_sqrt = sqrt(a * b);
        a = 0.5 * (a + b);
        b = ab_sqrt;

        d += d;
    }

    double temp = (atan(t) + mod * MY_PI) / (d * a);
    if (sign < 0) {
        temp = -temp;
    }
    return temp + npio2 * K;
}



// ---- ellpj: Jacobi elliptic functions sn, cn, dn, am ----
 GPUFUN 
void ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph)
{
    // Domain check
    if (m < 0.0 || m > 1.0) {
        *sn = *cn = *dn = *ph = 0.0;
        return;
    }

    // ---- Case 1: m extremely small (sinusoidal limit)
    if (m < 1.0e-9) {
        double t = sin(u);
        double c = cos(u);
        double ai = 0.25 * m * (u - t * c);

        *sn = t - ai * c;
        *cn = c + ai * t;
        *ph = u - ai;
        *dn = 1.0 - 0.5 * m * t * t;
        return;
    }

    // ---- Case 2: m extremely close to 1 (tanh limit)
    if (m > 0.9999999999) {
        double ai = 0.25 * (1.0 - m);
        double b = cosh(u);
        double t = tanh(u);
        double phi = 1.0 / b;
        double twon = b * sinh(u);

        *sn = t + ai * (twon - u) / (b * b);
        *ph = 2.0 * atan(exp(u)) - MY_PIO2 + ai * (twon - u) / b;

        ai *= t * phi;
        *cn = phi - ai * (twon - u);
        *dn = phi + ai * (twon + u);
        return;
    }

    // ---- General case: AGM iteration
    double a[9], c[9];
    a[0] = 1.0;
    c[0] = sqrt(m);
    double b = sqrt(1.0 - m);

    double twon = 1.0;
    int i = 0;

    while (fabs(c[i] / a[i]) > MY_MACHEP) {
        if (i >= 7) break;

        double ai = a[i];
        double bi = b;

        ++i;
        c[i] = 0.5 * (ai - bi);
        double ab = sqrt(ai * bi);
        a[i] = 0.5 * (ai + bi);
        b = ab;

        twon *= 2.0;
    }

    // phi = u * twon * a[i]
    double phi = twon * a[i] * u;

    // descending Landen transform
    while (i > 0) {
        double t = c[i] * sin(phi) / a[i];
        phi = 0.5 * (asin(t) + phi);
        --i;
    }

    double s = sin(phi);

    *sn = s;
    *cn = cos(phi);
    *dn = sqrt(1.0 - m * s * s);
    *ph = phi;
}


#endif // MYFUNCTIONS_H
