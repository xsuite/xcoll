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

// so this struct is not needed as the parameters will come from traj and seg. 
typedef struct {
    double s0;
    double x0;
    double theta;
    double s1;
    double x1;
    // double R;
    // double theta;
    double s2;
    double x2;
} Params;


double F_G(double t, Params* p){
    // return p->x1 + p->R*sin(t) - p->x0-(p->s1 - p->s0 + p->R*cos(t))*tan(p->theta);
    return (1-t)*p->x1+t*p->x2 - p->x0-((1-t)*p->s1+t*p->s2 - p->s0)*tan(p->theta);
}

double deriv_FG(double t, Params* p){
    // return p->R*cos(t) + p->R*sin(t)*tan(p->theta);
    return -p->x1+p->x2 + (p->s1-p->s2)*tan(p->theta);
}

double newton(double (*F_G)(double,Params*),double (*deriv_FG)(double,Params*), 
             double guess_t, Params* p) {
    for (int i = 0; i < 5000; i++) {
        double f = F_G(guess_t, p);
        double f_prime = deriv_FG(guess_t, p);
        if (fabs(f_prime) < 1e-10){
            return guess_t;
        } 
        double guess_new = guess_t - f / f_prime;
        if (fabs(guess_new - guess_t) < 1e-14){
            return guess_new;
        }
        guess_t = guess_new;
    }
    return guess_t;
}

void grid_search_and_newton(double (*F_G)(double,Params*), double (*deriv_FG)(double,Params*), double t_min, double t_max, double* roots, double max_crossings, int8_t* number_of_roots, Params* p) {
    /// Find the intervals where the function changes sign within the range [s_min, s_max]
    //  in which later Newton's method can be applied to find the root(s) for each interval
    double grid_step = (t_max - t_min) / 1000;
    int interval_count = 0;

    double prev_t   = t_min;
    double prev_val = F_G(prev_t, p);

    for (int i = 1; i <= 5000 - 1; i++) {
        if (interval_count >= max_crossings) break; // you cannot have more intervals than roots
        double curr_t = t_min + i * grid_step;
        double curr_val = F_G(curr_t, p);
        if (prev_val * curr_val < 0) {
            double initial_guess = 0.5 * (prev_t + curr_t);      // initial guess is midpoint
            roots[interval_count] = newton(F_G,deriv_FG, initial_guess, p);
            interval_count++;
        }
        prev_t = curr_t;
        prev_val = curr_val;
    }
    *number_of_roots = interval_count;
}

// // for testing
// int main(){
//     double s0 = 0.0;
//     double x0 = 1.5;
//     double s1 = 1.0;
//     double x1 = -1.0;
//     double s2 = 2.0;
//     double x2 = 2.0;
//     // double R = 3.0;
//     double theta = 20.0*3.14/180.0;
//     Params params = {s0, x0, theta, s1, x1, s2, x2};
//     double roots[1];
// //     // Define roots array and parameters
// //     Params_Line params = {Xo, *Ax, x2, s2, x1, s1};
// //     double roots[XC_LINE_CROSSINGS];
//     int8_t number_of_roots = 0;
//     grid_search_and_newton(F_G,deriv_FG, -3, 3, roots, 2, &number_of_roots, &params);

//     for (int i = 0; i < 1; ++i) {
//         if (roots[i] >= -3 && roots[i] <= 3) {
//             printf("The root t is %f\n", roots[i]);
//             printf("the other root g(t) is %f\n", ((1-roots[i])*s1+roots[i]*s2-s0)/cos(theta));
//             printf("xS: %f\n", (1-roots[i])*x1 + roots[i]*x2);
//             printf("xT: %f\n", x0 +sin(theta)*((1-roots[i])*s1+roots[i]*s2-s0)/cos(theta));
//             printf("sS: %f\n", (1-roots[i])*s1 + roots[i]*s2);
//             printf("sT: %f\n", s0 +((1-roots[i])*s1+roots[i]*s2-s0));
//         }
//     }
// }

#endif /* XCOLL_GEOM_FIND_ROOT_H */