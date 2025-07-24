// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_CHANNELLINGDEV_H
#define XCOLL_CHANNELLINGDEV_H


#include <math.h>
#include <stdio.h>

// Moliere functions 

static inline double E_T_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i, double a_TF) {
    return bpc*theta*theta / 2.0+U_N*(cosh(beta_i/a_TF*x) - 1.0);
}

static inline double m1_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i, double a_TF) {
    double E_T = E_T_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    return 1.0 + 2.0*U_N /E_T;
}

static inline double m1p_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i, double a_TF) {
    double E_T = E_T_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    return E_T/(E_T + 2.0*U_N);
}

static inline double nu_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i, double a_TF) {
    double E_T = E_T_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    return -copysign(1.0,theta)*beta_i/a_TF*sqrt(E_T/(2.0*bpc));
}

static inline double phi_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i, double a_TF) {
    double m1 = m1_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double m1p = m1p_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double arg_am = sqrt(m1) * tanh(x * beta_i / (2.0 * a_TF));
    if (arg_am > 1.0) arg_am = 1.0 - 1e-14;
    if (arg_am < -1.0) arg_am = -1.0 + 1e-14;
    double phi_in = asin(arg_am);
    double F = ellik(phi_in, m1p);
    return -sqrt(m1p) * F;
}

static inline double x_simplemoliere(double z, double x, double theta, double bpc, double U_N, double beta_i, double a_TF) {
    double m1 = m1_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double m1p = m1p_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double nu = nu_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double phi = phi_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double u = sqrt(m1)*(nu*z + phi);
    double sn, cn, dn, am;
    ellpj(u, m1p, &sn, &cn, &dn, &am);
    return -2.0 * a_TF / beta_i * asinh(sqrt(m1p) * sn / dn);
}

static inline double theta_simplemoliere(double z, double x, double theta, double bpc, double U_N, double beta_i, double a_TF) {
    double m1 = m1_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double m1p = m1p_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double nu = nu_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double phi = phi_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF);
    double u = sqrt(m1) * (nu * z + phi);
    double sn, cn, dn, am;
    ellpj(u, m1p, &sn, &cn, &dn, &am);
    return -2.0 * a_TF / beta_i * nu * cn / dn;
}
/*gpufun*/
void ChannellingDev_track_local_particle(ChannellingDevData el, LocalParticle* part0) {
    double length = ChannellingDevData_get_length(el);

    double beta_i = 0.573481;
    double a_TF = 0.194E-10;
    double U_N = 3.526347;
    double bpc = 150E9; 
    //start_per_particle_block (part0->part)
        double x0 = LocalParticle_get_x(part);
        double theta0 = LocalParticle_get_xp(part);
        double z = length;

        double x = x_simplemoliere(z, x0, theta0, bpc, U_N, beta_i, a_TF);
        double theta = theta_simplemoliere(z, x0, theta0, bpc, U_N, beta_i, a_TF);

        LocalParticle_set_x(part, x);
        LocalParticle_set_xp(part, theta);
    //end_per_particle_block
}


#endif /* XCOLL_CHANNELLINGDEV_H */
