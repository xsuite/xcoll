// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#define XC_NEWTON_EPSILON 1.e-10
#define XC_NEWTON_MAX_ITER 100  // Maximum number of iterations
#define XC_NEWTON_DERIVATIVE_TOL 1e-10  // Threshold for small derivative
// #define XC_GRID_MAX_INTER 10     // Maximum number of intervals for grid search
#define XC_GRID_POINTS 1000     // Number of points to search in grid

#include <math.h>
#include <stdio.h>

/*gpufun*/
// void grid_search_and_newton(LocalTrajectory traj, LocalSegment seg,
//                             double s_min, double s_max, double* roots, double max_crossings,
//                             int8_t* number_of_roots) {
//     /// Find the intervals where the function changes sign within the range [s_min, s_max]
//     //  in which later Newton's method can be applied to find the root(s) for each interval
//     double grid_step = (s_max - s_min) / XC_GRID_POINTS;
//     int interval_count = 0;

//     double prev_s   = s_min;
//     double prev_val = LocalTrajectory_func(traj, prev_s) - LocalSegment_func(seg, prev_s);

//     for (int i = 1; i <= XC_GRID_POINTS - 1; i++) {
//         if (interval_count >= max_crossings) break; // you cannot have more intervals than roots
//         double curr_s = s_min + i * grid_step;
//         double curr_val = LocalTrajectory_func(traj, curr_s) - LocalSegment_func(seg, curr_s);
//         if (prev_val * curr_val < 0) {
//             double initial_guess = 0.5 * (prev_s + curr_s);      // initial guess is midpoint
//             roots[interval_count] = newton(traj, seg, initial_guess);
//             interval_count++;
//         }
//         prev_s = curr_s;
//         prev_val = curr_val;
//     }
//     *number_of_roots = interval_count;
// }

// /*gpufun*/
// double newton(LocalTrajectory traj, LocalSegment seg, double guess_s) {
//     for (int i = 0; i < XC_NEWTON_MAX_ITER; i++) {
//         double f = LocalTrajectory_func(traj, guess_s) - LocalSegment_func(seg, guess_s);
//         double f_prime = LocalTrajectory_deriv(traj, guess_s) - LocalSegment_deriv(seg, guess_s);
//         if (fabs(f_prime) < XC_NEWTON_DERIVATIVE_TOL) return guess_s;

//         double guess_new = guess_s - f / f_prime;
//         if (fabs(guess_new - guess_s) < XC_NEWTON_EPSILON) return guess_new;
//         guess_s = guess_new;
//     }
//     return guess_s;
// }

// // --------------------------------------------------------------------------------------------

// Dont actually need this, just to make it easier to read. l and t need to be passed. 
typedef struct {
    double s0;
    double x0;
    double theta;
    double s1;
    double x1;
    double s2;
    double x2;
    double l;
    double t;
} Params;


void F_G(double FG[2], Params* p){
    // return p->x1 + p->R*sin(t) - p->x0-(p->s1 - p->s0 + p->R*cos(t))*tan(p->theta);
    printf("The t and l are,t = %f, l = %f\n", p->t, p->l);
    FG[0] = (1-p->t)*p->s1+p->t*p->s2 - (p->s0+p->l*cos(p->theta));
    FG[1] = (1-p->t)*p->x1+p->t*p->x2 - (p->x0+p->l*sin(p->theta));
    // Here we get the expr from segment and traj, and we connect them to create F1 and F2
    // F1 = sS(t) - sT(l) 
    // F2 = xS(t) - xT(l)


    // double sS = (1-p->t)*p->s1+p->t*p->s2;
    // double sT = (p->s0+p->l*cos(p->theta));
    // double xS =  (1-p->t)*p->x1+p->t*p->x2;
    // double xT = (p->x0+p->l*sin(p->theta));
    // printf("The functions are, sS = %f, sT = %f\n", sS, sT);
    // printf("The functions are, xS = %f, xT = %f\n", xS, xT);
    // printf("The functions are, f1 = %f, f2 = %f\n", FG[0], FG[1]);
    // return (1-t)*p->x1+t*p->x2 - p->x0-((1-t)*p->s1+t*p->s2 - p->s0)*tan(p->theta);
}

void get_inv_J(double J_inv[2][2], int8_t* no_crossing, Params* p){
    double J[2][2];
    // get derivatives
    J[0][0] = -p->s1 + p->s2; // get deriv. sS
    J[0][1] = -cos(p->theta); // get deriv. sT
    J[1][0] = -p->x1 + p->x2; // get deriv. xS
    J[1][1] = -sin(p->theta); // get deriv. xT
    
    // printf("Jacobian matrix:\n");
    // printf("[ %f, %f ]\n", J[0][0], J[0][1]);
    // printf("[ %f, %f ]\n", J[1][0], J[1][1]);
    double det = J[0][0] * J[1][1] - J[1][0]*J[0][1];
    if (det < XC_NEWTON_DERIVATIVE_TOL){
        printf("Determinant = 0! There is no crossing. \n");
        *no_crossing = 1;
        return;
    }
    J_inv[0][0] =  J[1][1] / det;
    J_inv[0][1] = -J[0][1] / det;
    J_inv[1][0] = -J[1][0] / det;
    J_inv[1][1] =  J[0][0] / det;
}

// needs rewrite but works
void newton(double guess_t, double guess_l, Params* p) {
    double J_inv[2][2];
    double FG[2];
    p->l = guess_l;
    p->t = guess_t;
    int8_t* no_crossing = 0;

    for (int i = 0; i < XC_NEWTON_MAX_ITER; i++) {
        F_G(FG, p);
        get_inv_J(J_inv, no_crossing, p); 
        if (no_crossing){

            return;
        }
        double new_t = p->t - (J_inv[0][0]*FG[0] + J_inv[0][1]*FG[1]);
        double new_l = p->l - (J_inv[1][0]*FG[0] + J_inv[1][1]*FG[1]);

    // Check for convergence
        if ((fabs(new_t -  p->t) < XC_NEWTON_DERIVATIVE_TOL) && (fabs(new_l - p->l) < XC_NEWTON_DERIVATIVE_TOL)){
            return;
        }
        // Update the guesses for the next iteration
        p->t = new_t;  // Keep p->t updated
        p->l = new_l;  // Keep p->l updated
    }
}

// int main() {
//     double s0 = 1.5;
//     double x0 = 0.0;
//     double theta = 20.0*M_PI/180.0;
//     double s1 = 1.0;
//     double x1 =-1.0;
//     double s2 =2.0;
//     double x2 =2.0;
//     double l = -1.0;
//     double t = -1.0;
//     Params p = {s0, x0, theta, s1, x1, s2, x2, l, t};
//     double guess_t = 0.5;
//     double guess_l = 0.5;
//     newton(&guess_t, &guess_l, &p);
//     return 0;
// }

#endif /* XCOLL_GEOM_FIND_ROOT_H */