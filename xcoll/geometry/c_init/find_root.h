// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#include <math.h>
#include <stdio.h>


typedef struct {
    double rC;     // length of position vector to first vertex
    double sin_tC; // angle of position vector to first vertex
    double cos_tC;
    double proj_l; // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    double proj_w; // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
    double l;      // length of the box
    double w;      // width of the box
    double sin_tb; // orientation of the box (angle of length wrt horizontal)
    double cos_tb;
} BoundingBox;


/*gpufun*/
int8_t bounding_boxes_overlap(BoundingBox* b1, BoundingBox* b2){
    // v1-v4 are the four vertices of the first box in counterclockwise order
    // w1-w4 are the four vertices of the second box in counterclockwise order
    // e1-e2 are the two axes of the first box
    // f1-f2 are the two axes of the second box
    double projs[4];
    double cos_tb1_tb2 = b1->cos_tb * b2->cos_tb + b1->sin_tb * b2->sin_tb;  // cos (tb1 - tb2)
    double sin_tb1_tb2 = b1->sin_tb * b2->cos_tb - b1->cos_tb * b2->sin_tb;  // sin (tb1 - tb2)

    // The length of the projection of vector W on E is given by: |W| cos(theta) where theta is the angle between W and E
    // projection of vertices of box 2 on the length axis of box 1 (e1)
    projs[0] = b2->rC * (b1->cos_tb * b2->cos_tC + b1->sin_tb * b2->sin_tC);  // first vertex w1:  |w1| cos (tb1 - tc2)
    projs[1] = projs[0] + b2->l * cos_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + b2->w * sin_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + b2->w * sin_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(b1->proj_l, b1->proj_l + b1->l, projs[0], projs[3])){ return false;}

    // length of projection of vertices of box 2 on the width axis of box 1 (e2)
    projs[0] = b2->rC * (b1->cos_tb * b2->sin_tC - b1->sin_tb * b2->cos_tC);  // first vertex w1:  |w1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] + b2->l * sin_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + b2->w * cos_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + b2->w * cos_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(b1->proj_w, b1->proj_w + b1->w, projs[0], projs[3])){ return false;}

    // length of projection of vertices of box 1 on the length axis of box 2 (f1)
    projs[0] = b2->rC * (b1->cos_tb * b2->cos_tC + b1->sin_tb * b2->sin_tC);  // first vertex w1:  |w1| cos (tb1 - tc2)
    projs[1] = projs[0] + b2->l * cos_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + b2->w * sin_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + b2->w * sin_tb1_tb2; // fourth vertex w4 = w1 + f2
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(b2->proj_l, b2->proj_l + b2->l, projs[0], projs[3])){ return false;}

    // length of projection of vertices of box 1 on the width axis of box 2 (f2)
    projs[0] = b2->rC * (b1->cos_tb * b2->sin_tC - b1->sin_tb * b2->cos_tC);  // first vertex w1:  |w1| cos(tb1 + pi/2 - tc2)
    projs[1] = projs[0] + b2->l * sin_tb1_tb2; // second vertex w2 = w1 + f1
    projs[2] = projs[1] + b2->w * cos_tb1_tb2; // third vertex w3 = w1 + f1 + f2
    projs[3] = projs[0] + b2->w * cos_tb1_tb2; // fourth vertex w4 = w1 + f2 
    sort_array_of_4_double(projs);
    if (!INTERVALS_OVERLAP(b2->proj_w, b2->proj_w + b2->w, projs[0], projs[3])){ return false;}

    return true;
}



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

void LocalCrossing_func(LocalTrajectory traj, LocalSegment seg, double TS[2], double l, double t){
    // Here we get the expr from segment and traj, and we connect them to create TSs and TSx
    // TSs = Ts(l) - Ss(t)
    // TSx = Tx(l) - Sx(t)
    TS[0] = LocalTrajectory_func_s(traj, l) - LocalSegment_func_s(seg, t);
    TS[1] = LocalTrajectory_func_x(traj, l) - LocalSegment_func_x(seg, t);
}

void LocalCrossing_inv_J(LocalTrajectory traj, LocalSegment seg, double J_inv[2][2], int8_t* no_crossing, double l, double t){
    double J[2][2];
    // get derivatives
    J[0][0] = LocalTrajectory_deriv_s(traj,l); // get deriv. dsT/dl LocalTrajectory_deriv_s
    J[0][1] = -LocalSegment_deriv_s(seg, t);     // get deriv. dsS/dt LocalSegment_deriv_s
    J[1][0] = LocalTrajectory_deriv_x(traj,l); // get deriv. dxT/dl LocalTrajectory_deriv_x
    J[1][1] = -LocalSegment_deriv_x(seg, t);     // get deriv. dxS/dt LocalSegment_deriv_x

    double det = J[0][0] * J[1][1] - J[1][0]*J[0][1];
    if (fabs(det) < XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL){
        printf("There is no crossing. \n");
        *no_crossing = 1;
        return;
    }
    J_inv[0][0] =  J[1][1] / det;
    J_inv[0][1] = -J[0][1] / det;
    J_inv[1][0] = -J[1][0] / det;
    J_inv[1][1] =  J[0][0] / det;
}

void newton(LocalTrajectory traj, LocalSegment seg, double* guess_l, double* guess_t, int8_t* no_crossing){
    double J_inv[2][2];
    double TS[2];

    for (int i = 0; i < XC_GEOM_ROOT_NEWTON_MAX_ITER; i++){
        F_G(traj, seg, TS, *guess_l, *guess_t);
        get_inv_J(traj, seg, J_inv, no_crossing, *guess_l, *guess_t);
        if (*no_crossing){
            return;
        }
        double new_t = *guess_t - (J_inv[0][0]*TS[0] + J_inv[0][1]*TS[1]);
        double new_l = *guess_l - (J_inv[1][0]*TS[0] + J_inv[1][1]*TS[1]);

        // Check for convergence
        if ((fabs(new_t -  *guess_t) < XC_GEOM_ROOT_NEWTON_EPSILON) && (fabs(new_l - *guess_l) < XC_GEOM_ROOT_NEWTON_EPSILON)){
            return;
        }
        // Update the guesses for the next iteration
        *guess_t = new_t;  // Keep *t updated
        *guess_l = new_l;  // Keep *l updated
    }
    *no_crossing = 1;
}

int8_t LocalCrossing_box_has_root(double TS_UL[2], double TS_UR[2], double TS_DL[2], double TS_DR[2]){

    // Evaluate F at corners
    evaluate_traj_seg(traj, seg, TS_DL, l1, t1);
    evaluate_traj_seg(traj, seg, TS_DR, l2, t1);
    evaluate_traj_seg(traj, seg, TS_UL, l1, t2);
    evaluate_traj_seg(traj, seg, TS_UR, l2, t2);

    // Check for sign changes for each component across the corners.
    // (For simplicity, here we check if the min and max differ in sign.)
    int8_t f1_has_sign_change = (min(TS_DL[0], TS_DR[0], TS_UL[0], TS_UR[0]) < 0 &&
                                max(TS_DL[0], TS_DR[0], TS_UL[0], TS_UR[0]) > 0);
    int8_t f2_has_sign_change = (min(TS_DL[1], TS_DR[1], TS_UL[1], TS_UR[1]) < 0 &&
                                max(TS_DL[1], TS_DR[1], TS_UL[1], TS_UR[1]) > 0);
    return (l1 < 0 && l2 > 0) || (l1 > 0 && l2 < 0) || (t1 < 0 && t2 > 0) || (t1 > 0 && t2 < 0);
}


void grid_search_and_newton(LocalTrajectory traj, LocalSegment seg, double* l, double* t, double l_min,
                            double l_max, double t_min, double t_max, double* roots_l, double* roots_t,
                            int max_crossings, int* number_of_roots){
    int N_l = 100;
    int N_t = 100;
    double grid_step_l = (l_max - l_min) / N_l; // this is just for now. We need to get interval for t and l.
    double grid_step_t = (t_max - t_min) / N_t;
    int n_roots = 0;
    double TS_prev, TS_curr;
    double prev_t = t_min;
    double prev_l = l_min;
    double curr_t;
    double curr_l;
    int8_t no_crossing = 0;

    LocalCrossing_func(traj, seg, TS_prev, prev_t, prev_l);

    // Stage 1: Coarse grid search over (l, t)
    for (int i=0; i < N_l-1; i++){
        double l1 = l_min + i * grid_step_l;
        double l2 = l_min + (i+1) * grid_step_l;
        for (int j=0; j < N_t-1; j++){
            double t1 = t_min + j * grid_step_t;
            double t2 = t_min + (j+1) * grid_step_t;
            F_G(traj, seg, FG_curr, curr_t, curr_l);





            if (n_roots >= max_crossings) {
                return;  // all possible roots have been found
            }
            curr_t = t_min + t_step;
            curr_l = l_min + l_step;
            F_G(traj, seg, FG_curr, curr_t, curr_l);

            if ((FG_prev[0] * FG_curr[0] < 0) || (FG_prev[1] * FG_curr[1] < 0)) {
                double initial_guess_t = 0.5*(curr_t - prev_t);
                double initial_guess_l = 0.5*(curr_l - prev_l);
                newton(traj, seg, &initial_guess_l, &initial_guess_t, &no_crossing);
                if (no_crossing){
                    printf("No crossing found with Newton's method. \n");
                } else{
                    roots_t[n_roots] = initial_guess_t;
                    roots_l[n_roots] = initial_guess_l;
                    n_roots++;
                }
            }
        }
    }
}

//  double prev_s   = s_min;
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