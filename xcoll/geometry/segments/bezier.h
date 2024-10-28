// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_BEZIERSEG_H
#define XCOLL_COLL_GEOM_BEZIERSEG_H

/*gpufun*/
void _hit_s_bezier(BezierSegment seg, double t, double multiplicity, int8_t* n_hit, double* s){
    double s1  = BezierSegment_get_s1(seg);
    double s2  = BezierSegment_get_s2(seg);
    double cs1 = BezierSegment_get_cs1(seg);
    double cs2 = BezierSegment_get_cs2(seg);
    double new_s = (1-t)*(1-t)*(1-t)*s1 + 3*(1-t)*(1-t)*t*cs1 + 3*(1-t)*t*t*cs2 + t*t*t*s2;
    for (int8_t i = 0; i < multiplicity; i++) {
        s[*n_hit] = new_s;
        (*n_hit)++;
    }
}

/*gpufun*/
void BezierSegment_crossing_drift(BezierSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
    // Get segment data
    double s1  = BezierSegment_get_s1(seg);
    double x1  = BezierSegment_get_x1(seg);
    double s2  = BezierSegment_get_s2(seg);
    double x2  = BezierSegment_get_x2(seg);
    double cs1 = BezierSegment_get_cs1(seg);
    double cx1 = BezierSegment_get_cx1(seg);
    double cs2 = BezierSegment_get_cs2(seg);
    double cx2 = BezierSegment_get_cx2(seg);
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
    if (fabs(a) < XC_EPSILON){
        if (fabs(b) < XC_EPSILON){
            if (fabs(c) < XC_EPSILON){
                if (fabs(d) < XC_EPSILON){
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
                    _hit_s_bezier(seg, t, 1, n_hit, s);
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
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
            }
        }
    } else {
        // Full cubic equation. Coefficients for the depressed cubic t^3 + p*t + q = 0:
        double p = (3*a*c - b*b)/(3*a*a);
        double q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
        double disc = -p*p*p/27 - q*q/4;  // This is the discriminant of the depressed cubic but divided by (4*27)
        if (fabs(disc) < XC_EPSILON){
            if (fabs(p) < XC_EPSILON){
                // One real root with multiplicity 3
                t = -b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 3, n_hit, s);
                }
            } else {
                // Two real roots (one simple and one with multiplicity 2)
                t = 3*q/p - b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
                t = -3*q/(2*p) - b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 2, n_hit, s);
                }
            }
        } else if (disc < 0){
            // One real root
            t = cbrt(-q/2 + sqrt(fabs(disc))) + cbrt(-q/2 - sqrt(fabs(disc))) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
        } else {
            // Three real roots
            double phi = acos(3*q/(2*p)*sqrt(fabs(3/p)));
            t = 2*sqrt(fabs(p/3))*cos(phi/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 2*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 4*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
        }
    }
}

#endif /* XCOLL_COLL_GEOM_BEZIERSEG_H */