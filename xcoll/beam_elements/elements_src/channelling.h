// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_CHANNELLINGDEV_H
#define XCOLL_CHANNELLINGDEV_H


#include <math.h>
#include <stdio.h>

// Moliere functions

/*gpufun*/
double U_simplemoliere(double x, double U_N, double beta_i, double a_TF) {
    return U_N*(cosh(beta_i/a_TF*x) - 1.0);
}

/*gpufun*/
double E_T_simplemoliere(double x, double theta, double bpc, double U_N,
                         double beta_i, double a_TF) {
    return bpc*theta*theta / 2.0 + U_simplemoliere(x, U_N, beta_i, a_TF);
}

/*gpufun*/
double m1_simplemoliere(double x, double theta, double bpc, double U_N,
                        double beta_i, double a_TF, double E_T ) {
    return 1.0 + 2.0*U_N /E_T;
}

/*gpufun*/
double m1p_simplemoliere(double x, double theta, double bpc, double U_N,
                         double beta_i, double a_TF, double E_T ) {
    return E_T/(E_T + 2.0*U_N);
}

/*gpufun*/
double nu_simplemoliere(double x, double theta, double bpc, double U_N,
                        double beta_i, double a_TF, double E_T) {
    // in unit 1.E4
    double sign = -1;
    if (theta <= 0) {
        sign = 1;
    }
    return sign*beta_i/a_TF*sqrt(E_T/(2.0*bpc));
}

/*gpufun*/
double phi_simplemoliere(double x, double theta, double nu, double bpc, double U_N,
                         double beta_i, double a_TF, double E_T) {
    double m = m1_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double mp = m1p_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double U = U_simplemoliere(x, U_N, beta_i, a_TF);
    double sign = -1;
    if (x <= 0) {
        sign = 1;
    }
    if (fabs(theta) < 1e-16) {
        // If theta is very close to 0, we can avoid numerical issues with asin(1)
        return sign*sqrt(mp)*ellpk(mp);
    }
    double phi_amplitude = asin(sqrt(m*U/(U+2*U_N)));
    // double alpha = 1/sqrt(m*pow(sinh(x*beta_i / (2*a_TF)), 2.) + pow(theta*beta_i /(2*a_TF*nu), 2.));
    // double phi_amplitude = asin(sqrt((1 - alpha*alpha)/mp));
    return sign*sqrt(mp)*ellik(phi_amplitude, mp);
}


/*gpufun*/
double x_simplemoliere(double z, double x, double theta, double nu, double bpc, double U_N,
                       double beta_i, double a_TF,double E_T, double phi) {
    // in Angstrom
    double m1 = m1_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double m1p = m1p_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double u = sqrt(m1) * (nu*1.E4 * z + phi);
    double sn, cn, dn, am;
    ellpj(u, m1p, &sn, &cn, &dn, &am);
    return -2.0 * a_TF / beta_i * asinh(sqrt(m1p) * sn / dn);
}

/*gpufun*/
double theta_simplemoliere(double z, double x, double theta, double nu, double bpc, double U_N,
                           double beta_i, double a_TF, double E_T, double phi) {
    // in urad
    double m1 = m1_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double m1p = m1p_simplemoliere(x, theta, bpc, U_N, beta_i, a_TF, E_T);
    double u = sqrt(m1) * (nu*1.E4 * z + phi);
    double sn, cn, dn, am;
    ellpj(u, m1p, &sn, &cn, &dn, &am);
    return -2.0 * a_TF / beta_i * nu * cn / dn;
}

/*gpufun*/
void ChannellingDev_track_local_particle(ChannellingDevData el, LocalParticle* part0) {
    double length = ChannellingDevData_get_length(el);

    double beta_i = 0.573481;
    double a_TF = 0.194;  // In Angstrom
    double U_N = 3.526347;
    double bpc = 0.150; // In TeV
    //start_per_particle_block (part0->part)
        double x0 = LocalParticle_get_x(part);
        double theta0 = LocalParticle_get_xp(part);
        double z = length;
        double E_T = E_T_simplemoliere(x0,  theta0,  bpc,  U_N,  beta_i,  a_TF);
        double nu = nu_simplemoliere(x0, theta0, bpc, U_N, beta_i, a_TF, E_T);
        double phi = phi_simplemoliere(x0,  theta0, nu,  bpc,  U_N,  beta_i,  a_TF, E_T);
        double x = x_simplemoliere(z, x0, theta0, nu, bpc, U_N, beta_i, a_TF, E_T, phi); // in Angstrom
        double theta = theta_simplemoliere(z, x0, theta0, nu, bpc, U_N, beta_i, a_TF, E_T, phi); // in urad

        LocalParticle_set_x(part, x);
        LocalParticle_set_xp(part, theta);
        // printf("ChannellingDev: x = %f, theta = %f, ET_before = %f, ET_after = %f\n", x*1.e10, theta*1.e9, E_T, E_T_simplemoliere(x,  theta,  bpc,  U_N,  beta_i,  a_TF));
    //end_per_particle_block
}


#endif /* XCOLL_CHANNELLINGDEV_H */
