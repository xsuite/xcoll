// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #


#ifndef BENT_CHANNELLING_KERNELS_H
#define BENT_CHANNELLING_KERNELS_H

#include <math.h>




// I am not sure which way is faster, I COMMENTED THEM OUT FOR NOW 

// beta_i / a_TF
//static inline double beta_over_aTF_(void) {
//    return beta_i/aTF;
//}

//static inline double aTF_over_beta_(void) {
//    return aTF/beta_i;
//}

//static inline double two_aTF_over_beta_(void) {
//    return 2.0*aTF/beta_i;
//}

// U_N = 2 Umax alpha_i / beta_i * exp( -beta_i dp / (2 aTF) )
static inline double U_N_(double Umax, double dp, 
    double alpha_i, double beta_i, double beta_over_aTF) {
    return 2.0 * Umax * alpha_i/beta_i
         * exp( -beta_over_aTF * (dp * 0.5) );
}


// Simplified Molière potential:
// U(x) = U_N (cosh(βx/aTF) − 1) = 2 U_N sinh²(βx / 2aTF)
static inline double U_simplemoliere(double x, double U_N, 
    double beta_over_aTF)
{   
    double t = beta_over_aTF*x*0.5;
    double sh = sinh(t);
    return 2.0*U_N*sh*sh;
}

// Total transverse energy E_T = p^2/(2βpc) + U(x)
static inline double E_T_simplemoliere(double theta, double bpc, double U)
{   
    double th2 = theta * theta;
    return 0.5*bpc*th2 + U; 
}




// m(E_T) = (E_T + 2U_N) / E_T
//static inline double m_simplemoliere(double E_T, double U_N)
//{
//    return 1.0 + 2.0 * U_N / E_T;
//}

// m' = 1/m
//static inline double mp_simplemoliere(double E_T, double U_N)
//{
//    return E_T / (E_T + 2.0 * U_N);
//}

// Nonlinear frequency parameter ν
static inline double nu_simplemoliere(double theta, double E_T, double bpc,
    double beta_over_aTF)
{
    double sign = (theta < 0.0) ? 1.0 : -1.0;

    return sign * beta_over_aTF * sqrt(E_T / (2.0 * bpc));
}

// Phase φ(x)
static inline double phi_simplemoliere(double x, double U_N, 
    double U, double m, double mp, double sqrt_mp)
{   
    double sign = (x < 0.0) ? 1.0 : -1.0;
    double denom = 1.0 / (U + 2.0 * U_N);
    double arg = sqrt(m * U * denom);
    // numerical safety
    if (arg > 1.0) arg = 1.0;
    if (arg < 0.0) arg = 0.0;
    double phi_ampl = asin(arg);

    return sign * sqrt_mp * ellik(phi_ampl, mp);
}



GPUFUN void motion_parameters_simplemoliere(
        double x0, double theta0, double z,
        double nu, double two_aTF_over_beta, double E_T,
        double phi,
        double m, double mp, double sqrt_m, double sqrt_mp,
        double* x, double* theta)
{
    if (E_T == 0.0) {
        *x     = x0;
        *theta = theta0;
        return;
    }

    double u = sqrt_m * (nu*z + phi);

    double sn, cn, dn, am;
    ellpj(u, mp, &sn, &cn, &dn, &am);
    // TO AVOID repeated division
    double inv_dn = 1.0/dn;
    *x     = -two_aTF_over_beta*asinh(sqrt_mp*sn*inv_dn);
    *theta = -two_aTF_over_beta*nu*cn*inv_dn;
}

#endif // BENT_CHANNELLING_KERNELS_H
