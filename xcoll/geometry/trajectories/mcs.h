// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_TRAJ_MCS_H
#define XCOLL_GEOM_TRAJ_MCS_H

#include <stdio.h>
#include <math.h>


#define MCS_AVERAGE_MOMENTUM 13.6e6
#define MCS_LOG_SCALE 3.821791440748616e-2
#define MCS_DERIV_LOG_SCALE 1.910895720374308e-2
#define MCS_DERIV_LOG_SHIFT 5.382179144074862e-1
#define MCS_RAN1_SCALE 2.886751345948129e-1  // 1/sqrt(12)
#define MCS_RAN2_SCALE 0.5


// /*gpufun*/
// void MultipleCoulombTrajectory_set_params(MultipleCoulombTrajectory traj, double X0,
//                                         double ran_1, double ran_2, LocalParticle part){
//     MultipleCoulombTrajectory_set_s0(traj, LocalParticle_get_s(part));
//     MultipleCoulombTrajectory_set_x0(traj, LocalParticle_get_s(part));
//     double xp = LocalParticle_get_exact_xp(part)
//     MultipleCoulombTrajectory_set_sin_t0(traj, xp / sqrt(1+xp*xp));
//     MultipleCoulombTrajectory_set_cos_t0(traj, 1. / sqrt(1+xp*xp));
//     MultipleCoulombTrajectory_set_tan_t0(traj, xp);
//     double beta = LocalParticle_get_rvv(part)*LocalParticle_get_beta0(part);
//     double q = LocalParticle_get_q0(part) * LocalParticle_get_charge_ratio(part);
//     MultipleCoulombTrajectory_set_Xt0(traj, X0*beta*beta / (q*q));
//     double pc = LocalParticle_get_p0c(part) * LocalParticle_get_charge_ratio(part) \
//                 / LocalParticle_get_chi(part) / LocalParticle_get_rpp(part); 
//     MultipleCoulombTrajectory_set_A0(traj, (ran_1*MCS_RAN1_SCALE + ran_2*MCS_RAN1_SCALE) * MCS_AVERAGE_MOMENTUM / pc);
//     MultipleCoulombTrajectory_set_B0(traj, ran_2 * MCS_AVERAGE_MOMENTUM / pc);
// }


/*gpufun*/
void MultipleCoulombTrajectory_set_params(MultipleCoulombTrajectory traj, double X0,
                                            double ran_1, double ran_2, double s0, double x0,
                                            double xp, double pc, double beta, double q){
    MultipleCoulombTrajectory_set_s0(traj, s0);
    MultipleCoulombTrajectory_set_x0(traj, x0);
    MultipleCoulombTrajectory_set_sin_t0(traj, xp / sqrt(1+xp*xp));
    MultipleCoulombTrajectory_set_cos_t0(traj, 1. / sqrt(1+xp*xp));
    MultipleCoulombTrajectory_set_tan_t0(traj, xp);
    MultipleCoulombTrajectory_set_Xt0(traj, X0*beta*beta / (q*q));
    MultipleCoulombTrajectory_set_A0(traj, (ran_1*MCS_RAN1_SCALE + ran_2*MCS_RAN2_SCALE) * MCS_AVERAGE_MOMENTUM / (beta*pc));
    MultipleCoulombTrajectory_set_B0(traj, ran_2 * MCS_AVERAGE_MOMENTUM / (beta*pc));
}

/*gpufun*/
double MultipleCoulombTrajectory_get_normalised_omega(MultipleCoulombTrajectory traj, double l){
    double Xt0 = MultipleCoulombTrajectory_get_Xt0(traj);  //  X0 𝛽^2 / q^2
    return sqrt(l/Xt0) * (1. + MCS_LOG_SCALE * log(l/Xt0));
}

/*gpufun*/
double MultipleCoulombTrajectory_get_normalised_omega_deriv(MultipleCoulombTrajectory traj, double l){
    double Xt0 = MultipleCoulombTrajectory_get_Xt0(traj);  //  X0 𝛽^2 / q^2
    return sqrt(l/Xt0) * (MCS_DERIV_LOG_SHIFT + MCS_DERIV_LOG_SCALE * log(l/Xt0)) / l;
}

/*gpufun*/
double MultipleCoulombTrajectory_func_s(MultipleCoulombTrajectory traj, double l){
    double s0 = MultipleCoulombTrajectory_get_s0(traj);
    double sin_t0 = MultipleCoulombTrajectory_get_sin_t0(traj);
    double cos_t0 = MultipleCoulombTrajectory_get_cos_t0(traj);
    double A0 = MultipleCoulombTrajectory_get_A0(traj);    // (𝜉1/√12 + 𝜉2/2) (13.6 MeV) / (pc)
    double omega_norm = MultipleCoulombTrajectory_get_normalised_omega(traj, l);
    return s0 + l*cos_t0 - l*A0*omega_norm*sin_t0;
}

/*gpufun*/
double MultipleCoulombTrajectory_func_x(MultipleCoulombTrajectory traj, double l){
    double x0 = MultipleCoulombTrajectory_get_x0(traj);
    double sin_t0 = MultipleCoulombTrajectory_get_sin_t0(traj);
    double cos_t0 = MultipleCoulombTrajectory_get_cos_t0(traj);
    double A0 = MultipleCoulombTrajectory_get_A0(traj);    // (𝜉1/√12 + 𝜉2/2) (13.6 MeV) / (pc)
    double omega_norm = MultipleCoulombTrajectory_get_normalised_omega(traj, l);
    return x0 + l*sin_t0 + l*A0*omega_norm*cos_t0;
}

/*gpufun*/
double MultipleCoulombTrajectory_func_xp(MultipleCoulombTrajectory traj, double l){
    double tan_t0 = MultipleCoulombTrajectory_get_tan_t0(traj);
    double B0 = MultipleCoulombTrajectory_get_B0(traj);    // 𝜉2 (13.6 MeV) / (pc)
    double omega_norm = MultipleCoulombTrajectory_get_normalised_omega(traj, l);
    return (tan_t0 + tan(B0*omega_norm)) / (1. - tan_t0*tan(B0*omega_norm));
}

/*gpufun*/
double MultipleCoulombTrajectory_deriv_s(MultipleCoulombTrajectory traj, double l){
    double sin_t0 = MultipleCoulombTrajectory_get_sin_t0(traj);
    double cos_t0 = MultipleCoulombTrajectory_get_cos_t0(traj);
    double A0 = MultipleCoulombTrajectory_get_A0(traj);    // (𝜉1/√12 + 𝜉2/2) (13.6 MeV) / (pc)
    double omega_norm = MultipleCoulombTrajectory_get_normalised_omega(traj, l);
    double omega_norm_deriv = MultipleCoulombTrajectory_get_normalised_omega_deriv(traj, l);
    return cos_t0 - A0*omega_norm*sin_t0 - l*A0*omega_norm_deriv*sin_t0;
}

/*gpufun*/
double MultipleCoulombTrajectory_deriv_x(MultipleCoulombTrajectory traj, double l){
    double sin_t0 = MultipleCoulombTrajectory_get_sin_t0(traj);
    double cos_t0 = MultipleCoulombTrajectory_get_cos_t0(traj);
    double A0 = MultipleCoulombTrajectory_get_A0(traj);    // (𝜉1/√12 + 𝜉2/2) (13.6 MeV) / (pc)
    double omega_norm = MultipleCoulombTrajectory_get_normalised_omega(traj, l);
    double omega_norm_deriv = MultipleCoulombTrajectory_get_normalised_omega_deriv(traj, l);
    return sin_t0 + A0*omega_norm*cos_t0 + l*A0*omega_norm_deriv*cos_t0;
}




// // MULTIPLE COULOMB SCATTERING VLIMIT ----------------------------------------------------------------------

// /*gpufun*/
// double MultipleCoulomb_y(double s, McsVlimitParams params){
//     double y  = McsVlimitParams_get_y(params);
//     double Ax = McsVlimitParams_get_Ax(params);
//     double Xo = McsVlimitParams_get_Xo(params);

//     double mcs = Ax * pow(sqrt(s/Xo),3.0) * (1.0/0.038 + log(s/Xo)); // just for clarity
//     return mcs - y;
// }

// /*gpufun*/
// double MultipleCoulombDeriv_y(double s, McsVlimitParams params){
//     double Ax = McsVlimitParams_get_Ax(params);
//     double Xo = McsVlimitParams_get_Xo(params);

//     double mcs_deriv = Ax/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo)); // just for clarity
//     return mcs_deriv;
// }
// /*gpufun*/
// int8_t MultipleCoulombTrajectory_vlimit(double* restrict_s, double s0, const double* Ax, const double Xo, double ymin, double ymax){
 
//     // int number_of_roots_min = 0;
//     // int number_of_roots_max = 0;
//     // double roots_min[1];
//     // double roots_max[1]; // can this be done smarter? 
//     // double s_max = 2.0;  // are we in this frame? I need coll length?  

//     // grid_search_and_newton(MultipleCoulomb_y, MultipleCoulombDeriv_y, s0, s_max, roots_max, 1, &params_max, &number_of_roots_max);
//     // grid_search_and_newton(MultipleCoulomb_y, MultipleCoulombDeriv_y, s0, s_max, roots_min, 1, &params_min, &number_of_roots_min);

//     // restrict_s[0] = roots_min[0];
//     // restrict_s[1] = roots_max[0];
//     // SWAP(restrict_s, 0, 1);   // To make sure these are sorted
//     // return 1;  // Default behavior: check overlap with horizontal crossings
//     // // }
// }

// // MULTIPLE COULOMB SCATTERING TRAJECTORY ------------------------------------------------------------------

// void generate_gaussian_random_numbers(double* z1, double* z2) {
//     // get z1 and z2 from Box-Muller transform
//     double v1, v2, r2;
//     do {
//         v1 = 2 * ((double)rand() / RAND_MAX) - 1;
//         v2 = 2 * ((double)rand() / RAND_MAX) - 1;
//         r2 = v1 * v1 + v2 * v2;
//     } while (r2 >= 1 || r2 == 0);
//     double a = sqrt(-2 * log(r2) / r2);
//     *z1 = v1 * a;
//     *z2 = v2 * a;
// }

// void A(double Xo, double p, double* Ax, double* Ay) {
//     // z1, z2 are gaussian random numbers. Xo is the radiation length, p is the momentum
//     double z1, z2;
//     generate_gaussian_random_numbers(&z1, &z2);
//     *Ax = (z1 / sqrt(12) + z2 / 2.0) * Xo * 13.6e-3 / p * 0.038;
//     *Ay = (z1 / sqrt(12) + z2 / 2.0) * Xo * 13.6e-3 / p * 0.038;
// }

// // void B(double* B) {
// //     *B = 1.0/0.038// + log(q**2/beta**2) // do we approx beta to one in the old code 
// // }

// // double MultipleCoulombTrajectory(double s, double const Xo, const double Ax){
// //     // MCS trajectory form PDG rewritted in terms of A, B and s/Xo. 
// //     return Ax * pow(sqrt(s/Xo),3.0) * (1.0/0.038 + log(s/Xo));
// // }

// // double MultipleCoulombTrajectory_prime_x(double s, double const Xo, const double Ax){ 
// //     return Ax/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo));
// // }

// // CURVE LENGTH ---------------------------------------------------------------------------------------------
// // Curve length is int_s1^s2 sqrt(1 + (dx/ds)^2 + (dy/ds)^2) ds
// double integrand_curve(double s, void *params) {
//     // Extract Ax, Ay, Xo from the params. Note: s is particle s.
//     double *p = (double *)params;
//     const double Ax = p[0];
//     const double Ay = p[1];
//     const double Xo = p[2];
//     double x;
//     double y;

//     x = Ax/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo));
//     y = Ay/Xo * (sqrt(s/Xo)*3.0/2.0*log(s/Xo)+1.0/0.038 + sqrt(s/Xo));
//     return sqrt(1 + x*x + y*y);
// }

// /*gpufun*/
// double MultipleCoulombTrajectory_length(double s0, double x0, double y0, const double s1, const double s2, 
//                                         const double* Ax, const double* Ay, const double Xo){
//     (void) s0;  // Avoid unused parameter warning
//     (void) x0;  // Avoid unused parameter warning
//     (void) y0;  // Avoid unused parameter warning, why do we have this though?
//     double n = 1000; // number of subintervals, must be even
//     double params[3] = {*Ax, *Ay, Xo}; // should i make Xo pointer for consistency?

//     // compare with adaptive gauss/kronrod quadrature later! 
//     double length = simpson(integrand_curve, s1, s2, n, params); // shold avoid magic numbers
//     return length;
// }

#endif /* XCOLL_GEOM_TRAJ_MCS_H */