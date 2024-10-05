// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_CROSSING_DRIFT_H
#define XCOLL_GEOM_CROSSING_DRIFT_H


// All segments have a function that calculates the s-coordinate of the crossing(s) with a
// particle trajectory in a drift (given by a line through (part_s, part_x) and a slope part_tan).
// The results are always stored in an array s, and n_hit keeps track of the number of hits over
// multiple segments.


// Line Segments
// -------------

/*gpufun*/
void _crossing_drift_line(int8_t* n_hit, double* s, double s1, double s2, double x1, double x2,
                          double part_s, double part_x, double part_tan){
    // We fill in the segment points in the particle trajectory equation; if the results have opposite sign,
    //the two points lie on different sides of the trajectory and hence the segment is crossed.
    double trajectory1 = x1 - part_x - (s1 - part_s)*part_tan;
    double trajectory2 = x2 - part_x - (s2 - part_s)*part_tan;
    if (trajectory1*trajectory2 <= 0){
        // It's a crossing
        if (fabs(s2 - s1) < XC_EPSILON){
            s[*n_hit] = s1;
            (*n_hit)++;
        } else {
            double poly_tan = (x2 - x1)/(s2 - s1);
            if (fabs(poly_tan - part_tan) < XC_EPSILON){
                // The trajectory is parallel to the segment.
                // TODO: this is situational; we should return s1 if get_s_first and current_s if after current_s
                s[*n_hit] = s1;
                (*n_hit)++;
            } else {
                // Normal crossing of two lines
                s[*n_hit] = (part_x - x1 + s1*poly_tan - part_s*part_tan)/(poly_tan - part_tan);
                (*n_hit)++;
            }
        }
    }
}

/*gpufun*/
void crossing_drift_line(void* segment, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    LineSegment seg = (LineSegment) segment;
    double s1 = seg->point1_s;
    double s2 = seg->point2_s;
    double x1 = seg->point1_x;
    double x2 = seg->point2_x;
    _crossing_drift_line(n_hit, s, s1, s2, x1, x2, part_s, part_x, part_tan);
}

/*gpufun*/
void crossing_drift_halfopenline(void* segment, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    HalfOpenLineSegment seg = (HalfOpenLineSegment) segment;
    double s1 = seg->point_s;
    double x1 = seg->point_x;
    // A half-open segment implies one of its points lies at +-inf.
    // In practice we just add a polygon point at the wall overflow (at 1000km for the x-coordinate).
    double x2 = 1.e6*seg->sign;
    double s2;
    if (fabs(seg->inv_slope) < XC_EPSILON){
        s2 = s1;
    } else {
        s2 = s1 - (x2 - x1)*seg->inv_slope;
    }
    _crossing_drift_line(n_hit, s, s1, s2, x1, x2, part_s, part_x, part_tan);
}


// Circular Segment
// ----------------

/*gpufun*/
void crossing_drift_circular(void* segment, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    CircularSegment seg = (CircularSegment) segment;
    double R = seg->R;
    double sC = seg->centre_s;
    double xC = seg->centre_x;
    double t1 = seg->point1_angle;
    double t2 = seg->point2_angle;
    // Calculate crossings
    int8_t reversed = 0;
    if (t2 < t1){
        reversed = 1;
    }
    double a = 1 + part_tan*part_tan;
    double bb = sC - part_tan*(part_x - xC - part_tan*part_s); // This is -b/2 with b from the quadratic formula
    double c = sC*sC + (part_x - xC - part_tan*part_s)*(part_x - xC - part_tan*part_s) - R*R;
    double disc = bb*bb - a*c; // This is  2*discriminant**2
    if (disc < 0){
        // No crossing
        return;
    }
    for (int8_t i = 0; i < 2; i++) {
        double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
        double new_s = (bb + sgnD*sqrt(fabs(disc))/a;
        double new_x = part_x + (new_s - part_s)*part_tan;
        double t = atan2(new_x - xC, new_s - sC);
        if (reversed){
            // t2 < t1, so we are looking at the inverted region of angles
            if (t1 >= t || t >= t2){
                s[*n_hit] = new_s;
                (*n_hit)++;
            }
        } else {
            if (t1 <= t && t <= t2){
                s[*n_hit] = new_s;
                (*n_hit)++;
            }
        }
    }
}


// Bézier Segment
// --------------

/*gpufun*/
double _hit_s_bezier(BezierSegment seg, double t, double multiplicity, int8_t* n_hit, double* s){
    double s1 = seg->point1_s;
    double cs1 = seg->control_point1_s;
    double s2 = seg->point2_s;
    double cs2 = seg->control_point2_s;
    double new_s = (1-t)*(1-t)*(1-t)*s1 + 3*(1-t)*(1-t)*t*cs1 + 3*(1-t)*t*t*cs2 + t*t*t*s2;
    for (int8_t i = 0; i < multiplicity; i++) {
        s[*n_hit] = new_s;
        (*n_hit)++;
    }
}

/*gpufun*/
void crossing_drift_bezier(void* segment, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    BezierSegment seg = (BezierSegment) segment;
    double s1 = seg->point1_s;
    double x1 = seg->point1_x;
    double cs1 = seg->control_point1_s;
    double cx1 = seg->control_point1_x;
    double s2 = seg->point2_s;
    double x2 = seg->point2_x;
    double cs2 = seg->control_point2_s;
    double cx2 = seg->control_point2_x;
    // The Bézier curve is defined by the parametric equations (with t in [0, 1]):
    // s(t) = (1-t)^3*s1 + 3(1-t)^2*t*cs1 + 3(1-t)*t^2*cs2 + t^3*s2
    // x(t) = (1-t)^3*x1 + 3(1-t)^2*t*cx1 + 3(1-t)*t^2*cx2 + t^3*x2
    double s0 = part_s;
    double x0 = part_x;
    double m = part_tan;
    // Plug the parametric eqs into the drift trajectory x(t) = m*(s(t) - s0) + x0 and solve for t
    // The solutions for t (which we get by Cardano's method) are valid if in [0, 1]
    double a = (m*s1 - x1) - (m*s2 - x2) - 3*(m*cs1 - cx1) + 3*(m*cs2 - cx2);
    double b = 6*(m*cs1 - cx1) - 3*(m*cs2 - cx2) - 3*(m*s1 - x1);
    double c = 3*(m*s1 - x1) - 3*(m*cs1 - cx1);
    double d = (m*s0 - x0) - (m*s1 - x1);
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


// Array of segments
// -----------------

/*gpufun*/
void crossing_drift(Segment* segments, int8_t n_segments, int8_t* n_hit, double* s, \
                    double part_s, double part_x, double part_tan){
    for (int8_t i=0; i<n_segments;i++) {
        int id = segments[i]->id;
        if (id == XC_LINESEGMENT_ID){
            crossing_drift_line(segments[i], n_hit, s, part_s, part_x, part_tan);
        } else if (id == XC_HALFOPENLINESEGMENT_ID){
            crossing_drift_halfopenline(segments[i], n_hit, s, part_s, part_x, part_tan);
        } else if (id == XC_CIRCULARSEGMENT_ID){
            crossing_drift_circular(segments[i], n_hit, s, part_s, part_x, part_tan);
        } else if (id == XC_BEZIERSEGMENT_ID){
            crossing_drift_bezier(segments[i], n_hit, s, part_s, part_x, part_tan);
        } // TODO: else throw fatal error
    }
    sort_array_of_double(s, (int64_t) *n_hit);
}

/*gpufun*/
void crossing_drift_vlimit(Segment* segments, int8_t n_segments, int8_t* n_hit, double* s, \
                           double part_s, double part_x, double part_tan_x, \
                           double part_y, double part_tan_y, \
                           double y_min, double y_max){
    if (fabs(part_tan_y) < XC_EPSILON){
        // Trajectory parallel to s axis
        if (part_y < y_min || part_y > y_max){
            // No crossing
            return;
        } else {
            // The particle is completely inside the vertical limits, so only check horizontal
            crossing_drift(segments, n_segments, n_hit, s, part_s, part_x, part_tan_x);
            return;
        }
    } else {
        crossing_drift(segments, n_segments, n_hit, s, part_s, part_x, part_tan_x);
        // restrict_s is the region [s0, s1] where the particle is inside the vertical limits
        double* restrict_s = (double*) malloc(2*sizeof(double));
        restrict_s[0] = (y_min - part_y)/part_tan_y + part_s;
        restrict_s[1] = (y_max - part_y)/part_tan_y + part_s;
        SWAP(restrict_s, 0, 1);   // To make sure these are sorted
        calculate_overlap_array_interval(s, n_hit, restrict_s);
        free(restrict_s);
    }
}

/*gpufun*/
int max_array_size_drift(Segment* segments, int8_t n_segments){
    int size = 0;
    for (int8_t i=0; i<n_segments;i++) {
        int id = segments[i]->id;
        if (id == XC_LINESEGMENT_ID){
            size += 1;
        } else if (id == XC_HALFOPENLINESEGMENT_ID){
            size += 1;
        } else if (id == XC_CIRCULARSEGMENT_ID){
            size += 2;
        } else if (id == XC_BEZIERSEGMENT_ID){
            size += 3;
        } // TODO: else throw fatal error
    }
    return size;
}

#endif /* XCOLL_GEOM_CROSSING_DRIFT_H */