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
#if 1
/* #ifdef __STDC__ */
#define ANSIPROT 1
int mtherr ( char *, int );
#else
int mtherr();
#endif

/* Variable for error reporting.  See mtherr.c.  */
extern int merror;


//--------------mtherr.c----------------

static int cephes_errno = 0;

/* Notice: the order of appearance of the following
 * messages is bound to the error codes defined
 * in mconf.h.
 */
static const char *ermsg[] = {
    "no error",
    "domain",       /* error code 1 */
    "singularity",  /* et seq.      */
    "overflow",
    "underflow",
    "total loss of precision",
    "partial loss of precision",
    "unknown"
};

/* @name is supposed to be the name of the function in
   which the error occurred; @code is an index into
   the array of error messages above; @arg is the offending
   argument, if @have_arg is non-zero, otherwise it is
   ignored.
*/

/*gpufun*/
int real_mtherr (char *name, int code, double arg,
			int have_arg)
{
    fprintf(stderr, "%s ", name);

    if (code <= 0 || code > PLOSS)
	code = PLOSS;

    cephes_errno = code;

    if (have_arg) {
	fprintf(stderr, "%s error (arg = %g)\n", ermsg[code], arg);
    } else {
	fprintf(stderr, "%s error\n", ermsg[code]);
    }

    return 0;
}

int mtherr_with_arg (char *name, int code, double arg)
{
    return real_mtherr(name, code, arg, 1);
}

int mtherr (char *name, int code)
{
    return real_mtherr(name, code, 0, 0);
}

int get_cephes_errno (void)
{
    int ret = cephes_errno;

    cephes_errno = 0; /* clear the code */

    return ret;
}







// IS THIS POLYNOMIAL EVALUATION ENOUGH? https://lost-contact.mit.edu/afs/lngs.infn.it/experiment/cuore/soft_node101/scipy-0.16.0/scipy/special/cephes/polevl.h
//---------------------------------------
// ---- Polynomial evaluation (Cephes) ----
/*gpufun*/
double polevl(double x, const double coef[], int N) {
    double ans = coef[0];
     for (int i = 1; i <= N; ++i) {
     ans = ans * x + coef[i];
    }
    return ans;
}

/*gpufun*/
double p1evl(double x, const double coef[], int N) {
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


/* UNK (portable double) – canonical values */
  const double ELLPK_P_UNK[11] = {
 1.37982864606273237150e-4, 2.28025724005875567385e-3,
 7.97404013220415179367e-3, 9.85821379021226008714e-3,
 6.87489687449949877925e-3, 6.18901033637687613229e-3,
 8.79078273952743772254e-3, 1.49380448916805252718e-2,
 3.08851465246711995998e-2, 9.65735902811690126535e-2,
 1.38629436111989062502e0
};
  const double ELLPK_Q_UNK[11] = {
 2.94078955048598507511e-5, 9.14184723865917226571e-4,
 5.94058303753167793257e-3, 1.54850516649762399335e-2,
 2.39089602715924892727e-2, 3.01204715227604046988e-2,
 3.73774314173823228969e-2, 4.88280347570998239232e-2,
 7.03124996963957469739e-2, 1.24999999999870820058e-1,
 4.99999999999999999821e-1
};

/* For completeness, expose DEC/IBMPC/MIEEE “variants” as doubles too.
   (They are numerically the same as UNK.) */
   double ELLPK_P_DEC[11]    = { /* = UNK */ 
 1.37982864606273237150e-4, 2.28025724005875567385e-3,
 7.97404013220415179367e-3, 9.85821379021226008714e-3,
 6.87489687449949877925e-3, 6.18901033637687613229e-3,
 8.79078273952743772254e-3, 1.49380448916805252718e-2,
 3.08851465246711995998e-2, 9.65735902811690126535e-2,
 1.38629436111989062502e0
};
   double ELLPK_Q_DEC[11]    = { /* = UNK */
 2.94078955048598507511e-5, 9.14184723865917226571e-4,
 5.94058303753167793257e-3, 1.54850516649762399335e-2,
 2.39089602715924892727e-2, 3.01204715227604046988e-2,
 3.73774314173823228969e-2, 4.88280347570998239232e-2,
 7.03124996963957469739e-2, 1.24999999999870820058e-1,
 4.99999999999999999821e-1
};

   double ELLPK_P_IBMPC[11]  = { /* = UNK */
 1.37982864606273237150e-4, 2.28025724005875567385e-3,
 7.97404013220415179367e-3, 9.85821379021226008714e-3,
 6.87489687449949877925e-3, 6.18901033637687613229e-3,
 8.79078273952743772254e-3, 1.49380448916805252718e-2,
 3.08851465246711995998e-2, 9.65735902811690126535e-2,
 1.38629436111989062502e0
};
  double ELLPK_Q_IBMPC[11]  = { /* = UNK */
 2.94078955048598507511e-5, 9.14184723865917226571e-4,
 5.94058303753167793257e-3, 1.54850516649762399335e-2,
 2.39089602715924892727e-2, 3.01204715227604046988e-2,
 3.73774314173823228969e-2, 4.88280347570998239232e-2,
 7.03124996963957469739e-2, 1.24999999999870820058e-1,
 4.99999999999999999821e-1
};

  double ELLPK_P_MIEEE[11]  = { /* = UNK */
 1.37982864606273237150e-4, 2.28025724005875567385e-3,
 7.97404013220415179367e-3, 9.85821379021226008714e-3,
 6.87489687449949877925e-3, 6.18901033637687613229e-3,
 8.79078273952743772254e-3, 1.49380448916805252718e-2,
 3.08851465246711995998e-2, 9.65735902811690126535e-2,
 1.38629436111989062502e0
};
   double ELLPK_Q_MIEEE[11]  = { /* = UNK */
 2.94078955048598507511e-5, 9.14184723865917226571e-4,
 5.94058303753167793257e-3, 1.54850516649762399335e-2,
 2.39089602715924892727e-2, 3.01204715227604046988e-2,
 3.73774314173823228969e-2, 4.88280347570998239232e-2,
 7.03124996963957469739e-2, 1.24999999999870820058e-1,
 4.99999999999999999821e-1
};

/* ---- Choose which set to use: define ELLPK_USE_{UNK,DEC,IBMPC,MIEEE} ---- */
#if defined(ELLPK_USE_DEC)
#  define ELLPK_P  ELLPK_P_DEC
#  define ELLPK_Q  ELLPK_Q_DEC
#elif defined(ELLPK_USE_IBMPC)
#  define ELLPK_P  ELLPK_P_IBMPC
#  define ELLPK_Q  ELLPK_Q_IBMPC
#elif defined(ELLPK_USE_MIEEE)
#  define ELLPK_P  ELLPK_P_MIEEE
#  define ELLPK_Q  ELLPK_Q_MIEEE
#else
#  define ELLPK_P  ELLPK_P_UNK     /* default */
#  define ELLPK_Q  ELLPK_Q_UNK
#endif

/* Optional: exact log(4) used in the small-x asymptotic path */
#define C1 1.3862943611198906188

double ellpk(x)
double x;
{

if( (x < 0.0) || (x > 1.0) )
	{
	mtherr( "ellpk", DOMAIN );
	return( 0.0 );
	}

if( x > MY_MACHEP )
	{
	return( polevl(x,ELLPK_P,10) - log(x) * polevl(x,ELLPK_Q,10) );
	}
else
	{
	if( x == 0.0 )
		{
		mtherr( "ellpk", SING );
		return( MY_MAXNUM );
		}
	else
		{
		return( C1 - 0.5 * log(x) );
		}
	}
}


// ---- ellik: incomplete elliptic integral of the first kind ----
/*gpufun*/
double ellik(double phi, double m) {
    if (m == 0.0)
        return phi;

    double a = 1.0 - m;
    if (a == 0.0) {
        if (fabs(phi) >= MY_PIO2)
            return MY_MAXNUM; 
        return log(tan((MY_PIO2 + phi) / 2.0));
    }

    int npio2 = (int)floor(phi / MY_PIO2);
    if (npio2 & 1) npio2 += 1;

    double K;
    if (npio2) {
        K = ellpk(a);
        phi = phi - npio2 * MY_PIO2;
    } else {
        K = 0.0;
    }

    int sign;
    if (phi < 0.0) {
        phi = -phi;
        sign = -1;
    } else {
        sign = 0;
    }

    double b = sqrt(a);
    double t = tan(phi);

    if (fabs(t) > 10.0) {
        // avoid deep recursion   
        double e = 1.0 / (b * t);
        if (fabs(e) < 10.0) {
            e = atan(e);
            if (npio2 == 0)
                K = ellpk(a);
            double temp = K - ellik(e, m);
            double out = (sign ? -1.0 : 1.0) * temp + npio2 * K;
            return out;
        }
    }

    a = 1.0;
    double c = sqrt(m);
    int d = 1, mod = 0;

    while (fabs(c / a) > MY_MACHEP) {
        double temp = b / a;
        phi = phi + atan(t * temp) + mod * MY_PI;
        mod = (int)((phi + MY_PIO2) / MY_PI);  

        t = t * (1.0 + temp) / (1.0 - temp * t * t);

        c = (a - b) / 2.0;

        double ab_sqrt = sqrt(a * b);
        //double a_old = a;
        a = (a + b) / 2.0;   
        b = ab_sqrt;         //  update b

        d += d;
    }

    double temp = (atan(t) + mod * MY_PI) / (d * a);
    if (sign < 0) temp = -temp;
    return temp + npio2 * K;
}


// ---- ellpj: Jacobi elliptic functions sn, cn, dn, am ----
 /*gpufun*/

void ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph) {
    if (m < 0.0 || m > 1.0) {
        *sn = *cn = *dn = *ph = 0.0;
    }
    // is it okay?
    if (m < 1.0e-9) {
        double t = sin(u), b = cos(u);
        double ai = 0.25 * m * (u - t * b);
        *sn = t - ai * b;
        *cn = b + ai * t;
        *ph = u - ai;
        *dn = 1.0 - 0.5 * m * t * t;
    }

    if (m >= 0.9999999999) {
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
    }

    double a[9], c[9], ai, b, phi, t, twon;
    a[0] = 1.0;
    b = sqrt(1.0 - m);
    c[0] = sqrt(m);
    twon = 1.0;
    int i = 0;

    while (fabs(c[i] / a[i]) > MY_MACHEP) {
        if (i > 7) break;                     // Cephes OVERFLOW 
        ai = a[i];
        ++i;
        c[i] = (ai - b) / 2.0;
        t = sqrt(ai * b);
        a[i] = (ai + b) / 2.0;
        b = t;
        twon *= 2.0;
    }

    phi = twon * a[i] * u;

    
    for (; i > 0; --i) {
        t = c[i] * sin(phi) / a[i];
        phi = 0.5 * (asin(t) + phi);
    }

    t = sin(phi);
    *sn = t;
    *cn = cos(phi);
    *dn = sqrt(1.0 - m * t * t);
    *ph = phi;
}

#endif // MYFUNCTIONS_H
