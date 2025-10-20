// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#include <math.h>
#include <stdio.h>

/*gpufun*/
void DriftTrajectory_analytical_solution_l(FindRoot finder, DriftTrajectory traj, double s, double x){
    if (DriftTrajectory_get_cos_t0(traj) - DriftTrajectory_get_cos_t0(traj) < XC_GEOM_EPSILON){
        FindRoot_set_solution_l(finder, FindRoot_get_num_solutions(finder), (x - DriftTrajectory_get_x0(traj)) / DriftTrajectory_get_sin_t0(traj));
    } else {
        FindRoot_set_solution_l(finder, FindRoot_get_num_solutions(finder), (s - DriftTrajectory_get_s0(traj)) / DriftTrajectory_get_cos_t0(traj));}
}

/*gpufun*/
void LineSegment_crossing_drift(FindRoot finder, LineSegment seg, DriftTrajectory traj, double s0, double x0, double xm){
    // Get segment data
    double s1 = LineSegment_get_s1(seg);
    double x1 = LineSegment_get_x1(seg);
    double s2 = LineSegment_get_s2(seg);
    double x2 = LineSegment_get_x2(seg);
    double denom = (x2 - x1) - (s2 - s1)*xm;
    if (fabs(denom) < XC_GEOM_EPSILON){
        // Trajectory is parallel to the segment
        if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_GEOM_EPSILON){
            // Trajectory overlaps with the segment
            // TODO: This is situational; we should return s1 if <get_s_first and current_s if after current_s
            //       For now we hit twice (because we go nor IN nor OUT)
            FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), 0.0);
            DriftTrajectory_analytical_solution_l(finder, traj, LineSegment_func_s(seg, 0.0), LineSegment_func_x(seg, 0.0));
            FindRoot_set_converged(finder, FindRoot_get_num_solutions(finder), 1);
            FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);

            FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), 1.0);
            DriftTrajectory_analytical_solution_l(finder, traj, LineSegment_func_s(seg, 1.0), LineSegment_func_x(seg, 1.0));
            FindRoot_set_converged(finder, FindRoot_get_num_solutions(finder), 1);
            FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
        } else {
            // No crossing
            return;
        }
    } else {
        double t = (x0 - x1 - (s0 - s1)*xm) / denom;
        if (t >= 0 && t <= 1){
            FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
            DriftTrajectory_analytical_solution_l(finder, traj, LineSegment_func_s(seg, t), LineSegment_func_x(seg, t));
            FindRoot_set_converged(finder, FindRoot_get_num_solutions(finder), 1);
            FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
        }
    }
}
/*gpufun*/
void HalfOpenLineSegment_crossing_drift(FindRoot finder, HalfOpenLineSegment seg, DriftTrajectory traj, double s0, double x0, double xm){
    // Get segment data
    double s1 = HalfOpenLineSegment_get_s1(seg);
    double x1 = HalfOpenLineSegment_get_x1(seg);
    double s2 = s1 + HalfOpenLineSegment_get_cos_t1(seg); //TODO: replace 10 with a parameter
    double x2 = x1 + HalfOpenLineSegment_get_sin_t1(seg);
    double denom = (x2 - x1) - (s2 - s1)*xm;
    if (fabs(denom) < XC_GEOM_EPSILON){
        // Trajectory is parallel to the segment
        if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_GEOM_EPSILON){
            FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), 0.0);
            DriftTrajectory_analytical_solution_l(finder, traj, HalfOpenLineSegment_func_s(seg, 0.0), HalfOpenLineSegment_func_x(seg, 0.0));
            FindRoot_set_converged(finder, FindRoot_get_num_solutions(finder), 1);
            FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
            // We do not add t=1 as it is a half-open segment
        } else {
            // No hit
            return;
        }
    } else {
        double t = (x0 - x1 - (s0 - s1)*xm) / denom;
        if (t >= 0){  // We do not check for t<=1 as it is a half-open segment
            FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
            DriftTrajectory_analytical_solution_l(finder, traj, HalfOpenLineSegment_func_s(seg, t), HalfOpenLineSegment_func_x(seg, t));
            FindRoot_set_converged(finder, FindRoot_get_num_solutions(finder), 1);
            FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
        }
    }
}
/*gpufun*/
void BezierSegment_crossing_drift(FindRoot finder, BezierSegment seg, DriftTrajectory traj, double s0, double x0, double xm){
    // Get segment data
    double s1  = BezierSegment_get__s1(seg);
    double x1  = BezierSegment_get__x1(seg);
    double s2  = BezierSegment_get__s2(seg);
    double x2  = BezierSegment_get__x2(seg);
    double cs1 = BezierSegment_get__cs1(seg);
    double cx1 = BezierSegment_get__cx1(seg);
    double cs2 = BezierSegment_get__cs2(seg);
    double cx2 = BezierSegment_get__cx2(seg);
    // The Bézier curve is defined by the parametric equations (with t in [0, 1]):
    // s(t) = (1-t)^3*s1 + 3(1-t)^2*t*cs1 + 3(1-t)*t^2*cs2 + t^3*s2
    // x(t) = (1-t)^3*x1 + 3(1-t)^2*t*cx1 + 3(1-t)*t^2*cx2 + t^3*x2
    // Plug the parametric eqs into the drift trajectory x(t) = m*(s(t) - s0) + x0 and solve for t
    // The solutions for t (which we get by Cardano's method) are valid if in [0, 1]
    double a = (xm*s1 - x1) - (xm*s2 - x2) - 3*(xm*cs1 - cx1) + 3*(xm*cs2 - cx2);
    double b = 6*(xm*cs1 - cx1) - 3*(xm*cs2 - cx2) - 3*(xm*s1 - x1);
    double c = 3*(xm*s1 - x1) - 3*(xm*cs1 - cx1);
    double d = (xm*s0 - x0) - (xm*s1 - x1);
    double t;
    // Edge cases
    if (fabs(a) < XC_GEOM_EPSILON){
        if (fabs(b) < XC_GEOM_EPSILON){
            if (fabs(c) < XC_GEOM_EPSILON){
                if (fabs(d) < XC_GEOM_EPSILON){
                    // The trajectory is on the Bézier curve
                    // TODO: This cannot happen because we don't have a cubic trajectory.
                    //       Adapt if these ever would be implemented.
                    return;
                } else {
                    // No solutions
                    return;
                }
            } else {
                // This is a linear equation
                t = -d/c;
                if (0 <= t && t <= 1){
                    FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                    DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                    FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
                }
            }
        } else {
            // This is a quadratic equation
            double disc = c*c - 4*b*d;
            if (disc < 0){
                // No solutions
                return;
            }
            for (int8_t i = 0; i < 2; i++) {
                double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
                t = (-c + sgnD*sqrt(fabs(disc)))/(2*b);
                if (0 <= t && t <= 1){
                    FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                    DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                    FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
                }
            }
        }
    } else {
        // Full cubic equation. Coefficients for the depressed cubic t^3 + p*t + q = 0:
        double p = (3*a*c - b*b)/(3*a*a);
        double q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
        double disc = -p*p*p/27 - q*q/4;  // This is the discriminant of the depressed cubic but divided by (4*27)
        if (fabs(disc) < XC_GEOM_EPSILON){
            if (fabs(p) < XC_GEOM_EPSILON){
                // One real root with multiplicity 3
                t = -b/(3*a);
                if (0 <= t && t <= 1){
                    FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                    DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                    FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
                }
            } else {
                // Two real roots (one simple and one with multiplicity 2)
                t = 3*q/p - b/(3*a);
                if (0 <= t && t <= 1){
                    FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                    DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                    FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
                }
                t = -3*q/(2*p) - b/(3*a);
                if (0 <= t && t <= 1){
                    FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                    DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                    FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
                }
            }
        } else if (disc < 0){
            // One real root
            t = cbrt(-q/2 + sqrt(fabs(disc))) + cbrt(-q/2 - sqrt(fabs(disc))) - b/(3*a);
            if (0 <= t && t <= 1){
                FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
            }
        } else {
            // Three real roots
            double phi = acos(3*q/(2*p)*sqrt(fabs(3/p)));
            t = 2*sqrt(fabs(p/3))*cos(phi/3) - b/(3*a);
            if (0 <= t && t <= 1){
                FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 2*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 4*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                FindRoot_set_solution_t(finder, FindRoot_get_num_solutions(finder), t);
                DriftTrajectory_analytical_solution_l(finder, traj, BezierSegment_func_s(seg, t), BezierSegment_func_x(seg, t));
                FindRoot_set_num_solutions(finder, FindRoot_get_num_solutions(finder)+1);
            }
        }
    }
}
/*gpufun*/
void LocalCrossing_func(LocalSegment seg, LocalTrajectory traj, double TS[2], double t, double l){
    // Here we get the expr from segment and traj, and we connect them to create TSs and TSx
    // TSs = Ts(l) - Ss(t)
    // TSx = Tx(l) - Sx(t)
    TS[0] = LocalTrajectory_func_s(traj, l) - LocalSegment_func_s(seg, t);
    TS[1] = LocalTrajectory_func_x(traj, l) - LocalSegment_func_x(seg, t);
}
/*gpufun*/
int8_t LocalCrossing_inv_J(LocalSegment seg, LocalTrajectory traj, double J_inv[2][2], double t, double l){
    double J[2][2];
    // get derivatives
    J[0][0] = LocalTrajectory_deriv_s(traj, l); // get deriv. dsT/dl LocalTrajectory_deriv_s
    J[0][1] = -LocalSegment_deriv_s(seg, t);     // get deriv. dsS/dt LocalSegment_deriv_s
    J[1][0] = LocalTrajectory_deriv_x(traj, l); // get deriv. dxT/dl LocalTrajectory_deriv_x
    J[1][1] = -LocalSegment_deriv_x(seg, t);     // get deriv. dxS/dt LocalSegment_deriv_x
    double det = J[0][0] * J[1][1] - J[1][0]*J[0][1];
    if (fabs(det) < XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL){
        return 0; // Failure, Jacobian is singular
    }
    J_inv[0][0] =  J[1][1] / det;
    J_inv[0][1] = -J[0][1] / det;
    J_inv[1][0] = -J[1][0] / det;
    J_inv[1][1] =  J[0][0] / det;
    return 1; // Success
}
/*gpufun*/
void FindRoot_newton(FindRoot finder, LocalSegment seg, LocalTrajectory traj, double guess_t, double guess_l,
                     int8_t num){
    // success = 0 means success, success = 1 means failure
    double J_inv[2][2];
    double TS[2];
    double corr0, corr1;
    double new_t, new_l;
    int8_t converged;

    for (int i = 0; i < XC_GEOM_ROOT_NEWTON_MAX_ITER; i++){
        LocalCrossing_func(seg, traj, TS, guess_t, guess_l);
        converged = LocalCrossing_inv_J(seg, traj, J_inv, guess_t, guess_l);
        if (converged != 1){
            printf("Jacobian is singular. Stopping Newton's method.\n");
            FindRoot_set_solution_t(finder, num, 1.e21);
            FindRoot_set_solution_l(finder, num, 1.e21);
            FindRoot_set_converged(finder, num, converged);
            continue;
        }
        // Residual convergence checks
        if (fabs(TS[0]) < XC_GEOM_ROOT_NEWTON_EPSILON && fabs(TS[1]) < XC_GEOM_ROOT_NEWTON_EPSILON) {
            FindRoot_set_solution_t(finder, num, guess_t);
            FindRoot_set_solution_l(finder, num, guess_l);
            FindRoot_set_converged(finder, num, converged);
            return;
        }
        corr0 = J_inv[0][0]*TS[0] + J_inv[0][1]*TS[1]; // delta for l
        corr1 = J_inv[1][0]*TS[0] + J_inv[1][1]*TS[1]; // delta for t
        new_t = guess_t - corr1;
        new_l = guess_l - corr0;
        // Check for parameter update convergence
        if ((fabs(new_t -  guess_t) < XC_GEOM_ROOT_NEWTON_EPSILON) && (fabs(new_l - guess_l) < XC_GEOM_ROOT_NEWTON_EPSILON)){
            FindRoot_set_solution_t(finder, num, guess_t);
            FindRoot_set_solution_l(finder, num, guess_l);
            FindRoot_set_converged(finder, num, converged);
            return;
        }
        // Update the guesses for the next iteration
        guess_t = new_t;  // Keep *t updated
        guess_l = new_l;  // Keep *l updated
    }
    printf("Warning: Newton's method did not converge within the maximum number of iterations (%d).\n", XC_GEOM_ROOT_NEWTON_MAX_ITER);
    printf("Last guess: t = %f, l = %f\n", guess_t, guess_l);
    FindRoot_set_solution_t(finder, num, 1.e21);
    FindRoot_set_solution_l(finder, num, 1.e21);
    FindRoot_set_converged(finder, num, converged);
}
// TODO: make this more clean pls
int8_t XC_SLICING_NUM_STEPS = 4;
int8_t XC_SLICING_MAX_NEST_LEVEL = 3;

/*gpufun*/
void slice_before_newton(FindRoot finder, LocalSegment seg, LocalTrajectory traj,
                         double t1, double t2, double l1, double l2, int8_t nest_level){
    // Prepare initial guesses for Newton-Raphson root finding
    BoundingBox_ _box_seg = {0};
    BoundingBox_ _box_traj = {0};
    BoundingBox box_seg  = &_box_seg;
    BoundingBox box_traj = &_box_traj;
    double t_step = (t2 - t1) / XC_SLICING_NUM_STEPS;
    double l_step = (l2 - l1) / XC_SLICING_NUM_STEPS;
    double t, l;
    int16_t num;
    for (int i = 0; i < XC_SLICING_NUM_STEPS; i++){
        t = t1 + i * t_step;
        LocalSegment_update_box(seg, box_seg, t, t + t_step);
        for (int j = 0; j < XC_SLICING_NUM_STEPS; j++){
            l = l1 + j * l_step;
            LocalTrajectory_update_box(traj, box_traj, l, l + l_step);
            if (BoundingBox_overlaps(box_seg, box_traj)){
                if (nest_level >= XC_SLICING_MAX_NEST_LEVEL - 1){
                    // We reached the maximum nesting level, return the midpoint of the current t-interval
                    if (FindRoot_get_num_solutions(finder) >= FindRoot_get_max_solutions(finder)){
                        printf("Warning: Maximum number of solutions (%d) reached in slice_before_newton. Some solutions may be missed.\n", FindRoot_get_max_solutions(finder));
                        fflush(stdout);
                        return;
                    }
                    num = FindRoot_get_num_solutions(finder);
                    printf("Guess found at i,j: %d, %d (nest level %d). t = %f, l = %f\n", i, j, nest_level, t + 0.5 * t_step, l + 0.5 * l_step);
                    fflush(stdout);
                    FindRoot_set_guess_t(finder, num, t + 0.5 * t_step);
                    FindRoot_set_guess_l(finder, num, l + 0.5 * l_step);
                    FindRoot_set_num_solutions(finder, num + 1);
                    t_step = (t2 - t) / XC_SLICING_NUM_STEPS; // reset t_step in order to avoid duplicates due to size of boxes
                    l_step = (l2 - l) / XC_SLICING_NUM_STEPS; // reset l_step
                } else {
                    slice_before_newton(finder, seg, traj, t, t + t_step, l, l + l_step, nest_level + 1);
                }
            }
        }
    }
    return;
}

/*gpufun*/
void find_crossing_approximate(FindRoot finder, LocalSegment seg, LocalTrajectory traj){
    printf("Finding crossing approximately...\n");
    slice_before_newton(finder, seg, traj, 0, 1, 0, 10, 0); // this does not work, its a bad solution beacause scaling changed other stuff
    // int8_t num_found = FindRoot_get_num_solutions(finder);
    for (int i = 0; i < FindRoot_get_num_solutions(finder); i++){
        double guess_t = FindRoot_get_guess_t(finder, i);
        double guess_l = FindRoot_get_guess_l(finder, i);
        FindRoot_newton(finder, seg, traj, guess_t, guess_l, i);
    }
}
/*gpufun*/
void FindRoot_find_crossing(FindRoot finder, LocalSegment seg, LocalTrajectory traj){ //, double guess_l, double guess_t){
// this will act as the crossing main function - C magic
// First we check for the analytical solutions by checking what trajectory we have
    printf("Finding crossing...\n");
    fflush(stdout);
    switch (LocalTrajectory_typeid(traj)){
        case LocalTrajectory_DriftTrajectory_t:
            printf("Using analytical crossing method for drift trajectory.\n");
            double sin = DriftTrajectory_get_sin_t0((DriftTrajectory) LocalTrajectory_member(traj));
            double cos = DriftTrajectory_get_cos_t0((DriftTrajectory) LocalTrajectory_member(traj));
            double xp  = sin / cos;
            double s0  = DriftTrajectory_get_s0((DriftTrajectory) LocalTrajectory_member(traj));
            double x0  = DriftTrajectory_get_x0((DriftTrajectory) LocalTrajectory_member(traj));
            switch (LocalSegment_typeid(seg)){
                case LocalSegment_LineSegment_t:
                    LineSegment_crossing_drift(finder, (LineSegment) LocalSegment_member(seg), (DriftTrajectory) LocalTrajectory_member(traj), s0, x0, xp);
                    return;
                    break;
                case LocalSegment_HalfOpenLineSegment_t:
                    return HalfOpenLineSegment_crossing_drift(finder, (HalfOpenLineSegment) LocalSegment_member(seg), (DriftTrajectory) LocalTrajectory_member(traj), s0, x0, xp);
                    break;
                case LocalSegment_BezierSegment_t:
                    return BezierSegment_crossing_drift(finder, (BezierSegment) LocalSegment_member(seg), (DriftTrajectory) LocalTrajectory_member(traj), s0, x0, xp);
                default:
                    // Custom segment
                    return find_crossing_approximate(finder, seg, traj);
            }
            break;
        default:
            printf("Using approximate crossing method.\n");
            return find_crossing_approximate(finder, seg, traj);
    }
}
#endif /* XCOLL_GEOM_FIND_ROOT_H */