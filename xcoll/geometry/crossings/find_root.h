// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#include <math.h>
#include <stdio.h>

/*gpufun*/
void LocalCrossing_func(LocalTrajectory traj, LocalSegment seg, double TS[2], double l, double t){
    // Here we get the expr from segment and traj, and we connect them to create TSs and TSx
    // TSs = Ts(l) - Ss(t)
    // TSx = Tx(l) - Sx(t)
    TS[0] = LocalTrajectory_func_s(traj, l) - LocalSegment_func_s(seg, t);
    TS[1] = LocalTrajectory_func_x(traj, l) - LocalSegment_func_x(seg, t);
}
/*gpufun*/
void LocalCrossing_inv_J(LocalTrajectory traj, LocalSegment seg, double J_inv[2][2], double success[1], double l, double t){
    double J[2][2];
    // get derivatives
    J[0][0] = LocalTrajectory_deriv_s(traj, l); // get deriv. dsT/dl LocalTrajectory_deriv_s
    J[0][1] = -LocalSegment_deriv_s(seg, t);     // get deriv. dsS/dt LocalSegment_deriv_s
    J[1][0] = LocalTrajectory_deriv_x(traj, l); // get deriv. dxT/dl LocalTrajectory_deriv_x
    J[1][1] = -LocalSegment_deriv_x(seg, t);     // get deriv. dxS/dt LocalSegment_deriv_x
    double det = J[0][0] * J[1][1] - J[1][0]*J[0][1];
    if (fabs(det) < XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL){
        success[0] = 1;
        return; 
    }
    J_inv[0][0] =  J[1][1] / det;
    J_inv[0][1] = -J[0][1] / det;
    J_inv[1][0] = -J[1][0] / det;
    J_inv[1][1] =  J[0][0] / det;
    }

/*gpufun*/
void find_root_newton(find_root finder, LocalTrajectory traj, LocalSegment seg, double guess_l, double guess_t){
    // success = 0 means success, success = 1 means failure
    double J_inv[2][2];
    double TS[2];
    double success[1] = {0};
    double corr0, corr1;
    double new_t, new_l;
    for (int i = 0; i < XC_GEOM_ROOT_NEWTON_MAX_ITER; i++){
        LocalCrossing_func(traj, seg, TS, guess_l, guess_t);
        LocalCrossing_inv_J(traj, seg, J_inv, success, guess_l, guess_t);
        if (success[0]){
            printf("Jacobian is singular. Stopping Newton's method.\n");
            return;
        }
        // Residual convergence check
        if (fabs(TS[0]) < XC_GEOM_ROOT_NEWTON_EPSILON && fabs(TS[1]) < XC_GEOM_ROOT_NEWTON_EPSILON) {
            find_root_set_solution_l(finder, guess_l);
            find_root_set_solution_t(finder, guess_t);
            return;
        }
        corr0 = J_inv[0][0]*TS[0] + J_inv[0][1]*TS[1]; // delta for l
        corr1 = J_inv[1][0]*TS[0] + J_inv[1][1]*TS[1]; // delta for t
        new_t = guess_t - corr1;
        new_l = guess_l - corr0;
        // Check for parameter update convergence
        if ((fabs(new_t -  guess_t) < XC_GEOM_ROOT_NEWTON_EPSILON) && (fabs(new_l - guess_l) < XC_GEOM_ROOT_NEWTON_EPSILON)){
            find_root_set_solution_l(finder, guess_l);
            find_root_set_solution_t(finder, guess_t); // TODO: remove this, and the find root class member later.
            return;
        }
        // Update the guesses for the next iteration
        guess_t = new_t;  // Keep *t updated
        guess_l = new_l;  // Keep *l updated
    }
    printf("Warning: Newton's method did not converge within the maximum number of iterations (%d).\n", XC_GEOM_ROOT_NEWTON_MAX_ITER);
    success[0] = 1;
}

void slicing_boxes(find_root finder, LocalTrajectory traj, LocalSegment seg){
    // We have MCS, and are now finding what segment we are dealing with.
    // cast by (LineSegment) LocalSegment_get_member(seg)
    double tol    = 0.1;
    double guess_t, guess_l;
    switch (LocalSegment_typeid(seg)){
        case LocalSegment_LineSegment_t:
            guess_l = MultipleCoulombTrajectory_prepare_newton((MultipleCoulombTrajectory) LocalTrajectory_member(traj), 
                                                                       LineSegment_getp_box((LineSegment) LocalSegment_member(seg)), tol);
            guess_t = LineSegment_prepare_newton((LineSegment) LocalSegment_member(seg),
                                                         MultipleCoulombTrajectory_getp_box((MultipleCoulombTrajectory) LocalTrajectory_member(traj)), tol);
            find_root_newton(finder, traj, seg, guess_l, guess_t);
            return;
            break;
        case LocalSegment_HalfOpenLineSegment_t:
            guess_l = MultipleCoulombTrajectory_prepare_newton((MultipleCoulombTrajectory) LocalTrajectory_member(traj), 
                                                                       HalfOpenLineSegment_getp_box((HalfOpenLineSegment) LocalSegment_member(seg)), tol);
            guess_t = HalfOpenLineSegment_prepare_newton((HalfOpenLineSegment) LocalSegment_member(seg), 
                                                                 MultipleCoulombTrajectory_getp_box((MultipleCoulombTrajectory) LocalTrajectory_member(traj)), tol);
            find_root_newton(finder, traj, seg, guess_l, guess_t);
            return;
            break;
        case LocalSegment_BezierSegment_t:
            guess_l = MultipleCoulombTrajectory_prepare_newton((MultipleCoulombTrajectory) LocalTrajectory_member(traj), 
                                                                       BezierSegment_getp_box((BezierSegment) LocalSegment_member(seg)), tol);
            guess_t = BezierSegment_prepare_newton((BezierSegment) LocalSegment_member(seg), 
                                                                 MultipleCoulombTrajectory_getp_box((MultipleCoulombTrajectory) LocalTrajectory_member(traj)), tol);
            find_root_newton(finder, traj, seg, guess_l, guess_t);
            return;
            break;
  // // case LocalSegment_CircularSegment_t:
  // //     return prepare_newton_circular(traj, seg);
  // //     break;
        default:
            // Custom segment
            return;
    }
}

void find_root_find_crossing(find_root finder, LocalTrajectory traj, LocalSegment seg){ //, double guess_l, double guess_t){
// this will act as the crossing main function - C magic
// First we check for the analytical solutions by checking what trajectory we have
    switch (LocalTrajectory_typeid(traj)){
        // case LocalTrajectory_DriftTrajectory_t:
        //     switch (LocalSegment_typeid(seg)){
        //         // may have to add #ifndef LOCALSEGMENT_SKIP_LINESEGMENT
        //         case LocalSegment_LineSegment_t:
        //             return crossing_drift_line(traj, seg);
        //             break;
        //         case LocalSegment_HalfOpenLineSegment_t:
        //             return crossing_drift_halfline(traj, seg);
        //             break;
        //         // case LocalSegment_CircularSegment_t:
        //         //     return crossing_drift_circular(traj, seg);
        //         //     break;
        //         case LocalSegment_BezierSegment_t:
        //             return crossing_drift_bezier(traj, seg);
        //             break;
        //         default:
        //             // Custom segment
        //             return crossing_boxes_newton(traj, seg);
        //     }
        //     break;
        // case LocalTrajectory_CircularTrajectory_t:
        //     switch (LocalSegment_typeid(seg)){
        //         case LocalSegment_LineSegment_t:
        //             return crossing_circ_line(traj, seg);
        //             break;
        //         case LocalSegment_HalfOpenLineSegment_t:
        //             return crossing_circ_halfline(traj, seg);
        //             break;
        //         // case LocalSegment_CircularSegment_t:
        //         //     return crossing_circ_circular(traj, seg);
        //         //     break;
        //         case LocalSegment_BezierSegment_t:
        //             return crossing_circ_bezier(traj, seg);
        //             break;
        //         default:
        //             // Custom segment
        //             return crossing_boxes_newton(traj, seg);
        //     }
        //     break;
    default:
        // No analytical code for MultipleCoulombTrajectory, we now move on to figure out what segment we have
        return slicing_boxes(finder, traj, seg); //, guess_l, guess_t);
    }
} 
#endif /* XCOLL_GEOM_FIND_ROOT_H */