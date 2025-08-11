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

static inline double m1_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i, double a_TF,double E_T ) {
    return 1.0 + 2.0*U_N /E_T;
}

static inline double m1p_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i, double a_TF,double E_T ) {
    return E_T/(E_T + 2.0*U_N);
}

static inline double nu_simplemoliere(double x, double theta, double bpc, double U_N, double beta_i, double a_TF,double E_T) {
    return -copysign(1.0,theta)*beta_i/a_TF*sqrt(E_T/(2.0*bpc));
}

static inline double phi_simplemoliere(double x, double theta,
                                       double bpc, double U_N,
                                       double beta_i, double a_TF,double E_T) {
    double arg = - (beta_i * x) / (2.0 * a_TF);
    double phi_amplitude = atan(sinh(arg));
    double m = E_T / (E_T + 2.0 * U_N);

    return ellik(phi_amplitude, m);
}


static inline double x_simplemoliere(double z, double x, double theta, double bpc, double U_N, double beta_i, double a_TF,double E_T) {
    double m1 = m1_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double m1p = m1p_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double nu = nu_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double phi = phi_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double u = sqrt(m1)*(nu*z + phi);
    double sn, cn, dn, am;
    ellpj(u, m1p, &sn, &cn, &dn, &am);
    return -2.0 * a_TF / beta_i * asinh(sqrt(m1p) * sn / dn);
}

static inline double theta_simplemoliere(double z, double x, double theta, double bpc, double U_N, double beta_i, double a_TF,double E_T) {
    double m1 = m1_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double m1p = m1p_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double nu = nu_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double phi = phi_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
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
        double E_T=E_T_simplemoliere(x0,  theta0,  bpc,  U_N,  beta_i,  a_TF);
        double x = x_simplemoliere(z, x0, theta0, bpc, U_N, beta_i, a_TF, E_T);
        double theta = theta_simplemoliere(z, x0, theta0, bpc, U_N, beta_i, a_TF, E_T);

        LocalParticle_set_x(part, x);
        LocalParticle_set_xp(part, theta);
    //end_per_particle_block
}


#endif /* XCOLL_CHANNELLINGDEV_H */
