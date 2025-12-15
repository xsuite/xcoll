// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef BENT_CHANNELLING_INTEGRATORS_H
#define BENT_CHANNELLING_INTEGRATORS_H

#include <stdint.h>  // for int8_t

// ============================================================================
// Yoshida coefficients 
// ============================================================================
static const double yoshida_chi1[] = {
    0.0, 0.0, 0.0, 0.0,
    1.351207191959658,   // order 4
    0.0,
    1.174671758089364,   // order 6
    0.0,
    1.116182939325386,   // order 8
    0.0,
    1.087027106299171,   // order 10
    0.0,
    1.069565719632538    // order 12
};

static const double yoshida_chi0[] = {
    0.0, 0.0, 0.0, 0.0,
    -1.702414383919315,  // order 4
    0.0,
    -1.349343516178727,  // order 6
    0.0,
    -1.232365878650771,  // order 8
    0.0,
    -1.174054212598341,  // order 10
    0.0,
    -1.139131439265076   // order 12
};

// ============================================================================
// BASE SIGNATURE (M2/M3 H0/H1)
// ============================================================================
#define FM_SIG \
    double length, double x0, double px0, \
    double bpc,  double U_N, double U, double sqrt_U, double R, \
    double* x, double* px

    #define FM_SIG_H1 \
    double length, double x0, double px0, \
    double bpc,  double U_N, double U, \
    double aTF_over_beta, double beta_over_aTF, double R,\
    double* x, double* px
// ============================================================================
// M4 EXTENDED SIGNATURE (Simplified Molière model)
// ============================================================================
#define FM4_SIG \
    double length, double x0, double px0, \
    double bpc, double U_N,  double beta_over_aTF, double two_aTF_over_beta,\
    double* x, double* px

#define FM4_SIG_H1 \
    double length, double x0, double px0, double R,\
    double* x, double* px

// ============================================================================
// BASE SIGNATURE (O(2) integrators of M2/M3
// ============================================================================
#define FM_SIG_INT \
    double length, double x0, double px0, \
    double bpc,  double U_N, double U, double sqrt_U, \
    double aTF_over_beta, double beta_over_aTF, double R,\
    double* x, double* px

// ============================================================================
// BASE SIGNATURE (O(2) integrator of M4
// ============================================================================


#define FM4_SIG_INT \
    double length, double x0, double px0, \
    double bpc, double U_N,  double beta_over_aTF, double two_aTF_over_beta,\
    double R,\
    double* x, double* px




// ============================================================================
// ------------------------ METHOD 2: H0 / H1 ------------------------
// ============================================================================
GPUFUN void fM2H0(FM_SIG){
    // Some variables are directly computed just before applying yoshida  
   
    double u  = sqrt_U * length;
    double cosu = cos(u);
    double sinu = sin(u);

    *x  = x0*cosu + px0/sqrt_U*sinu;
    *px = -x0*sqrt_U*sinu + px0*cosu;
}
GPUFUN void fM2H1(FM_SIG_H1){
      
    *x  = x0;                 
    *px = px0 + U*x0*length - U*aTF_over_beta*sinh(beta_over_aTF*x0)*length - length/R;    // kick
}

// ============================================================================
// ------------------------ METHOD 3: H0 / H1 ------------------------
// ============================================================================
GPUFUN void fM3H0(FM_SIG) {
      
    double u = sqrt_U * length;

    double inv_UR  = 1.0 / ( (sqrt_U*sqrt_U) * R );
    double inv_sUR = 1.0 / ( sqrt_U * R );

    double cosu = cos(u);
    double sinu = sin(u);

    *x  = (x0 + inv_UR)*cosu + (px0/sqrt_U)*sinu - inv_UR;
    *px = (-x0*sqrt_U - inv_sUR)*sinu + px0*cosu;
}

GPUFUN void fM3H1(FM_SIG_H1) {
      
    *x = x0;
    *px = px0 + (U*x0 - U*aTF_over_beta*sinh(beta_over_aTF*x0))*length;

}
// ============================================================================
// ------------------------ METHOD 4: H0 / H1 ------------------------
// ============================================================================

GPUFUN void fM4H0(FM4_SIG)
{
    double U = U_simplemoliere(x0, U_N, beta_over_aTF);
    double E_T = E_T_simplemoliere(px0, bpc, U);

    if (E_T <= 0.0) { *x=x0; *px=px0; return; }

    double m = (E_T + 2.0 * U_N) / E_T;
    double mp      = 1.0/m;
    double sqrt_m = sqrt(m);
    double sqrt_mp = 1/sqrt_m;
    double nu = nu_simplemoliere(px0, E_T, bpc, beta_over_aTF);
    double phi = phi_simplemoliere(x0, U_N, U, m, mp, sqrt_mp);

    double x_new, th_new;
// Still needs to be optimised
    motion_parameters_simplemoliere(
         x0, px0, length,
         nu, two_aTF_over_beta, E_T, phi,
         m, mp, sqrt_m, sqrt_mp,
         &x_new, &th_new);

    *x  = x_new;
    *px = th_new;
}


GPUFUN void fM4H1(FM4_SIG_H1) {
    // Depends only on R

    if (R == 0.0) {*x  = x0;*px = px0;
        return;
    }

    *x  = x0;
    *px = px0 - length/R;
}

// ============================================================================
// ------------------------ M2 O(2) Integrators ------------------------
// ============================================================================
GPUFUN void fM2O2v1(FM_SIG_INT) {
    double half = 0.5 * length;
    double x1, px1, x2, px2;

    fM2H0(half, x0, px0, bpc, U_N, U, sqrt_U, R, &x1, &px1);
    fM2H1(length, x1, px1, bpc, U_N, U, aTF_over_beta, beta_over_aTF,
        R, &x2, &px2);
    fM2H0(half, x2, px2, bpc, U_N, U, sqrt_U, R, x, px);
}

GPUFUN void fM2O2v2(FM_SIG_INT) {
    double half = 0.5 * length;
    double x1, px1, x2, px2;

    fM2H1(half, x0, px0, bpc, U_N, U, aTF_over_beta, beta_over_aTF,
        R, &x1, &px1);
    fM2H0(length, x1, px1, bpc, U_N, U, sqrt_U, R, &x2, &px2);
    fM2H1(half, x2, px2,  bpc, U_N, U, aTF_over_beta, beta_over_aTF,
        R, x, px);
}

// ============================================================================
// ------------------------ M3 O(2) Integrators ------------------------
// ============================================================================
GPUFUN void fM3O2v1(FM_SIG_INT) {
    double half = 0.5 * length;
    double x1, px1, x2, px2;

    fM3H0(half, x0, px0, bpc, U_N, U, sqrt_U, R, &x1, &px1);
    fM3H1(length, x1, px1, bpc, U_N, U, aTF_over_beta, beta_over_aTF, R,
        &x2, &px2);
    fM3H0(half, x2, px2, bpc, U_N, U, sqrt_U, R, x,  px);
}

GPUFUN void fM3O2v2(FM_SIG_INT) {
    double half = 0.5 * length;
    double x1, px1, x2, px2;

    fM3H1(half, x0, px0, bpc, U_N, U, aTF_over_beta, beta_over_aTF, R,
        &x1, &px1);
    fM3H0(length, x1, px1, bpc, U_N, U, sqrt_U, R, &x2, &px2);
    fM3H1(half, x2, px2,bpc, U_N, U, aTF_over_beta, beta_over_aTF, R,
        x, px);
}

// ============================================================================
// ------------------------ M4 O(2) Integrators ------------------------
// ============================================================================

GPUFUN void fM4O2v1(FM4_SIG_INT) {
    double half = 0.5 * length;
    double x1, px1, x2, px2;

    // H0(h/2)
    fM4H0(half, x0, px0, bpc, U_N, beta_over_aTF, two_aTF_over_beta,
        &x1, &px1);
    // H1(h)
    fM4H1(length, x1, px1, R, &x2, &px2);
    // H0(h/2)
    fM4H0(half, x2, px2, bpc, U_N, beta_over_aTF, two_aTF_over_beta,
         x, px);
}

GPUFUN void fM4O2v2(FM4_SIG_INT) {
    double half = 0.5 * length;
    double x1, px1, x2, px2;

    // H1(h/2)
    fM4H1(
        half, x0, px0, R, &x1, &px1);
    // H0(h)
    fM4H0(length, x1, px1, bpc, U_N, beta_over_aTF, two_aTF_over_beta,
        &x2, &px2);
    // H1(h/2)
    fM4H1(half, x2, px2, R, x, px);
}

// ============================================================================
// ------------------------ YOSHIDA WRAPPERS (M2/M3/M4) ------------------------
// ============================================================================

GPUFUN void fM2_apply_yoshida(
    double length, double Umax, double aTF, double alpha_i, double beta_i,
    double dp, double x0, double px0, double bpc, double R,
    int8_t order, int8_t variant,
    double* x, double* px)
{    
    double aTF_over_beta = aTF/beta_i;
    double beta_over_aTF = 1/aTF_over_beta;
    double U_N = U_N_(Umax, dp, alpha_i, beta_i, beta_over_aTF);
    double U  = (U_N / bpc)*beta_over_aTF*beta_over_aTF;  
    double sqrt_U = sqrt(U);

    if (order <= 2) {
        (variant == 1 ? fM2O2v1 : fM2O2v2)(
            length, x0, px0, bpc, U_N, U,sqrt_U,aTF_over_beta, beta_over_aTF,
             R, x, px);
        return;
    }

    double chi1 = yoshida_chi1[order];
    double chi0 = yoshida_chi0[order];
    double x1, px1, x2, px2;

    fM2_apply_yoshida(
        chi1 * length, Umax, aTF, alpha_i, beta_i, dp,
        x0, px0, bpc, R, order - 2, variant, &x1, &px1
    );

    fM2_apply_yoshida(
        chi0 * length, Umax, aTF, alpha_i, beta_i, dp,
        x1, px1, bpc, R, order - 2, variant, &x2, &px2
    );

    fM2_apply_yoshida(
        chi1 * length, Umax, aTF, alpha_i, beta_i, dp,
        x2, px2, bpc, R, order - 2, variant, x, px
    );
}

GPUFUN void fM3_apply_yoshida(
    double length, double Umax, double aTF, double alpha_i, double beta_i,
    double dp, double x0, double px0, double bpc, double R,
    int8_t order, int8_t variant,
    double* x, double* px)
{   
    double aTF_over_beta = aTF/beta_i;
    double beta_over_aTF = 1/aTF_over_beta;
    double U_N = U_N_(Umax, dp, alpha_i, beta_i, beta_over_aTF);
    double U  = (U_N / bpc)*beta_over_aTF*beta_over_aTF;  
    double sqrt_U = sqrt(U);

    if (order <= 2) {
        (variant == 1 ? fM3O2v1 : fM3O2v2)(
            length, x0, px0, bpc, U_N, U,sqrt_U,aTF_over_beta, beta_over_aTF,
              R, x, px);
        return;
    }

    double chi1 = yoshida_chi1[order];
    double chi0 = yoshida_chi0[order];
    double x1, px1, x2, px2;

    fM3_apply_yoshida(
        chi1 * length, Umax, aTF, alpha_i, beta_i, dp,
        x0, px0, bpc, R, order - 2, variant, &x1, &px1
    );

    fM3_apply_yoshida(
        chi0 * length,  Umax, aTF, alpha_i, beta_i, dp,
        x1, px1, bpc, R, order - 2, variant, &x2, &px2
    );

    fM3_apply_yoshida(
        chi1 * length,  Umax, aTF, alpha_i, beta_i, dp,
        x2, px2, bpc, R, order - 2, variant, x, px
    );
}

GPUFUN void fM4_apply_yoshida(
    double length, double Umax, double aTF, double alpha_i, double beta_i,
    double dp,double x0, double px0,
    double bpc, double R, int8_t order, int8_t variant,
    double* x, double* px)
{   
    double beta_over_aTF = beta_i/aTF;
    double two_aTF_over_beta = 2/beta_over_aTF;
    double U_N = U_N_(Umax, dp, alpha_i, beta_i, beta_over_aTF);

    if (order <= 2) {
        if (variant == 1) {
            fM4O2v1(
                length, x0, px0, bpc, U_N, beta_over_aTF,two_aTF_over_beta,
                R, x, px);
        }
        else {
            fM4O2v2(
                length, x0, px0,bpc, U_N, beta_over_aTF, two_aTF_over_beta,
                R, x, px);
        }
        return;
    }

    // ------------------------------
    double chi1 = yoshida_chi1[order];
    double chi0 = yoshida_chi0[order];

    double s1 = chi1 * length;
    double s0 = chi0 * length;

    double x1, px1;
    double x2, px2;

    
    fM4_apply_yoshida(s1,  Umax, aTF, alpha_i, beta_i, dp,
        x0, px0, bpc, R, order - 2, variant,
        &x1, &px1
    );

    fM4_apply_yoshida(
        s0,  Umax, aTF, alpha_i, beta_i, dp,
        x1, px1, bpc, R, order - 2, variant,
        &x2, &px2
    );

    fM4_apply_yoshida(
        s1,  Umax, aTF, alpha_i, beta_i, dp,
        x2, px2, bpc, R, order - 2, variant,
        x, px
    );
}

#endif // BENT_CHANNELLING_INTEGRATORS_H
