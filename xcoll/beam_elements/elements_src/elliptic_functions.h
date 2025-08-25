#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

#include <math.h>

// ---- Mathematical constants (renamed to avoid conflicts) ----
#define MY_PI       3.14159265358979323846
#define MY_PIO2     1.57079632679489661923
//Cephes uses MACHEP 2^(-53) as a tight bound on relative rounding error in many formulas and stopping criteria.
#define MY_MACHEP   1.11022302462515654042E-16
#define MY_MAXNUM   1.7976931348623158E308


// IS THIS POLYNOMIAL EVALUATION ENOUGH???
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
/*gpufun*/
double ellpk(double x) {
    // log(4) constant used in the small-x asymptotic
    const double MY_LOG4 = 1.3862943611198906188; // = log(4)

    // Domain check, returns 0.0 on domain error)
    if (x < 0.0 || x > 1.0) {
        return 0.0;
    }

    // Small-x handling !!!!!!
    if (x <= MY_MACHEP) {
        if (x == 0.0) {
            // Singular at m -> 1 (x -> 0)
            return MY_MAXNUM;
        } else {
            // Asymptotic for very small x
            return MY_LOG4 - 0.5 * log(x);
        }
    }

    // Polynomial approximation: P(x) - log(x) Q(x)
    static const double P[] = {
        1.37982864606273237150E-4,
        2.28025724005875567385E-3,
        7.97404013220415179367E-3,
        9.85821379021226008714E-3,
        6.87489687449949877925E-3,
        6.18901033637687613229E-3,
        8.79078273952743772254E-3,
        1.49380448916805252718E-2,
        3.08851465246711995998E-2,
        9.65735902811690126535E-2,
        1.38629436111989062502E0
    };
    static const double Q[] = {
        2.94078955048598507511E-5,
        9.14184723865917226571E-4,
        5.94058303753167793257E-3,
        1.54850516649762399335E-2,
        2.39089602715924892727E-2,
        3.01204715227604046988E-2,
        3.73774314173823228969E-2,
        4.88280347570998239232E-2,
        7.03124996963957469739E-2,
        1.24999999999870820058E-1,
        4.99999999999999999821E-1
    };

    if (x == 1.0) {
        // K(0) = pi/2  ( m = 0)
        return MY_PIO2;
    }

    return polevl(x, P, 10) - log(x) * polevl(x, Q, 10);
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
