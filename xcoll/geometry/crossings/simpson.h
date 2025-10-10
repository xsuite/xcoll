// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SIMPSON_H
#define XCOLL_GEOM_SIMPSON_H
#include <stdio.h>
#include <math.h>
// Simpson's Rule function for numerical integration
// Simpsons rule approximates the integral of a function by dividing the integration interval into an even number of subintervals
// and applying a weighted sum of function values at the endpoints and midpoints of these subintervals.
// It returns the approximate value of the integral.
// (b - a) / (3 * subintervals) * (f(l1) + 4f(subinterval1) + 2f(even subinterval) + 4f(odd subinterval)) + ... + f(l2))
void simpson(FindRoot finder, LocalTrajectory traj, int subintervals) {
    double l1 = 0;
    double l2 = FindRoot_get_solution_l(finder, 0);
    if (subintervals % 2 != 1) {
        subintervals++;                                         // Requires an even number of subintervals
    }
    //double l1 = LocalTrajectory_get_l1(traj); l1 is always 0
    double step_size = (l2 - l1) / subintervals;
    double sum = sqrt(1 + LocalTrajectory_deriv_x(traj, l1)*LocalTrajectory_deriv_x(traj, l1)) + 
                 sqrt(1 + LocalTrajectory_deriv_x(traj, l2)*LocalTrajectory_deriv_x(traj, l2));  // f(l1) + f(l2)
    // Add subintervals to the sum
    for (int i = 1; i < subintervals; i++) {
        double l = l1 + i * step_size;
        // Even indices (except the endpoints) are multiplied by 2. Odd indices are multiplied by 4
        if (i % 2 == 0) {
            sum += 2.0 * sqrt(1 + LocalTrajectory_deriv_x(traj, l)*LocalTrajectory_deriv_x(traj, l));
        } else {
            sum += 4.0 * sqrt(1 + LocalTrajectory_deriv_x(traj, l)*LocalTrajectory_deriv_x(traj, l));
        }
    }
    sum *= step_size / 3;
    FindRoot_set_path_length(finder, sum);
    return;
}
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
#endif /* XCOLL_GEOM_SIMPSON_H */