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
double U_simplemoliere(double x, double U_N, double beta_i_over_a_TF) {
    //return U_N*(cosh(beta_i_over_a_TF*x) - 1.0);
    // In order to be stable near zero, I use the other form
    double sinh_;
    sinh_=sinh(beta_i_over_a_TF*x/2);
    return 2*U_N*sinh_*sinh_;
}

/*gpufun*/
double E_T_simplemoliere(double x, double theta, double bpc, double U_N,
                         double beta_i_over_a_TF) {
    return bpc*theta*theta / 2.0 + U_simplemoliere(x, U_N, beta_i_over_a_TF);
}

/*gpufun*/
double m_simplemoliere(double U_N, double E_T ) {
    return 1.0 + 2.0*U_N /E_T;
}

/*gpufun*/
double mp_simplemoliere(double U_N, double E_T ) {
    return E_T/(E_T + 2.0*U_N);
}

/*gpufun*/
double nu_simplemoliere(double theta, double bpc, double beta_i_over_a_TF, double E_T) {
    // in unit 1.E4
    double sign = -1;
    if (theta < 0) {
        sign = 1;
        }
    return sign*beta_i_over_a_TF*sqrt(E_T/(2.0*bpc));
}

/*gpufun*/
double phi_simplemoliere(double x, double theta, double nu, double bpc, double U_N,
                         double beta_i_over_a_TF, double m, double mp, double sqrt_mp) {
    double U = U_simplemoliere(x, U_N, beta_i_over_a_TF);
    double sign = -1;
    if (x < 0) {
        sign = 1;
    }
    //if (fabs(theta)< 1e-14) {
        // If theta is very close to 0, we can avoid numerical issues with asin(1)
        //return sign*sqrt_mp*ellpk(mp);
    //}
    double phi_amplitude = asin(sqrt(m*U/(U+2*U_N)));
    // double alpha = 1/sqrt(m*pow(sinh(x*beta_i / (2*a_TF)), 2.) + pow(theta*beta_i /(2*a_TF*nu), 2.));
    // double phi_amplitude = asin(sqrt((1 - alpha*alpha)/mp));
    return sign*sqrt_mp*ellik(phi_amplitude, mp);
}

/*gpufun*/
void motion_parameters(double x0, double theta0, double z, double nu, double E_T, double U_N, double twotimes_a_TF_over_beta_i, double phi, double m, double mp, double sqrt_mp,
    /*out*/ double* x_out, /*out*/ double* theta_out) {
    if (E_T == 0) {
        *x_out = 0;
        *theta_out=0;
    }
    else {
        //double u = sqrt(m) * (nu * 1.0e4 * z + phi);
        double u = sqrt(m)*(nu*z + phi);
        double sn, cn, dn, am;
        ellpj(u, mp, &sn, &cn, &dn, &am);
    
        *x_out = -twotimes_a_TF_over_beta_i * asinh(sqrt_mp * sn / dn);
        *theta_out = -twotimes_a_TF_over_beta_i * nu * cn / dn;
    }
    //if (fabs(theta0) > 1e-1) {
     //   *x_out = -twotimes_a_TF_over_beta_i * asinh(sqrt_mp * sn / dn);
   // } 
    //else {
    //    *x_out=twotimes_a_TF_over_beta_i/2*acosh(1+E_T/U_N);
    //}
}




// /*gpufun*/
//double x_simplemoliere(double z, double x, double theta, double nu,
 //                      double beta_i, double a_TF, double phi, double m, double mp) {
    // in Angstrom
    //double u = sqrt(m) * (nu*1.E4 * z + phi);
    //double sn, cn, dn, am;
    //ellpj(u, mp, &sn, &cn, &dn, &am);
    //return -2.0 * a_TF / beta_i * asinh(sqrt(mp) * sn / dn);
//}

///*gpufun*/
//double theta_simplemoliere(double z, double x, double theta, double nu,
 //                          double beta_i, double a_TF, double phi, double m, double mp) {
    // in urad
   // double u = sqrt(m) * (nu*1.E4 * z + phi);
   // double sn, cn, dn, am;
   // ellpj(u, mp, &sn, &cn, &dn, &am);
  //  return -2.0 * a_TF / beta_i * nu * cn / dn;
//}

/*gpufun*/
void ChannellingDev_track_local_particle(ChannellingDevData el, LocalParticle* part0) {
    double length = ChannellingDevData_get_length(el);

    double beta_i = 0.573481;

    //double a_TF = 0.194;  // In Angstrom
    // !!!!!
    double a_TF = 0.194e-10;
    // more constants, avoid calculations over and over 
    double twotimes_a_TF_over_beta_i = 2.0*a_TF/beta_i;
    double beta_i_over_a_TF = beta_i/a_TF;

    double U_N = 3.526347;
    // !!!!!
    //double bpc = 0.150; // In TeV
    double bpc = 150e9;
    //start_per_particle_block (part0->part)
        double x0 = LocalParticle_get_x(part);
        double theta0 = LocalParticle_get_xp(part);
        double z = length;
        double E_T = E_T_simplemoliere(x0,  theta0,  bpc,  U_N, beta_i_over_a_TF);
        double m = m_simplemoliere(U_N, E_T);
        double mp = mp_simplemoliere(U_N, E_T);

        //new constant
        double sqrt_mp = sqrt(mp);


        double nu = nu_simplemoliere(theta0, bpc, beta_i_over_a_TF, E_T);
        double phi = phi_simplemoliere(x0,  theta0, nu,  bpc,  U_N,  beta_i_over_a_TF, m, mp, sqrt_mp);
        double x, theta;
        motion_parameters(x0, theta0, z, nu, E_T, U_N, twotimes_a_TF_over_beta_i, phi, m, mp, sqrt_mp, &x, &theta);
        //double x = x_simplemoliere(z, x0, theta0, nu, beta_i, a_TF, phi, m, mp); // in Angstrom
        //double theta = theta_simplemoliere(z, x0, theta0, nu, beta_i, a_TF, phi, m, mp); // in urad


    



        LocalParticle_set_x(part, x);
        LocalParticle_set_xp(part, theta);
        //printf("ChannellingDev: x = %f, theta = %f, nu = %f, ET_before = %f", x, theta, nu, E_T);
    //end_per_particle_block
}


#endif /* XCOLL_CHANNELLINGDEV_H */
