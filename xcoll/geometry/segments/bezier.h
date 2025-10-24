// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_TRAJ_BEZIER_H
#define XCOLL_GEOM_TRAJ_BEZIER_H


/*gpufun*/
int8_t BezierSegment_bounded_below(BezierSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
int8_t BezierSegment_bounded_above(BezierSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
double BezierSegment_func_s(BezierSegment seg, double t){
    double s1  = BezierSegment_get__s1(seg);
    double s2  = BezierSegment_get__s2(seg);
    double cs1 = BezierSegment_get__cs1(seg);
    double cs2 = BezierSegment_get__cs2(seg);
    return (1-t)*(1-t)*(1-t)*s1 + 3*(1-t)*(1-t)*t*cs1 + 3*(1-t)*t*t*cs2 + t*t*t*s2;
}

/*gpufun*/
double BezierSegment_func_x(BezierSegment seg, double t){
    double x1  = BezierSegment_get__x1(seg);
    double x2  = BezierSegment_get__x2(seg);
    double cx1 = BezierSegment_get__cx1(seg);
    double cx2 = BezierSegment_get__cx2(seg);
    return (1-t)*(1-t)*(1-t)*x1 + 3*(1-t)*(1-t)*t*cx1 + 3*(1-t)*t*t*cx2 + t*t*t*x2;
}

/*gpufun*/
double BezierSegment_deriv_s(BezierSegment seg, double t){
    double s1  = BezierSegment_get__s1(seg);
    double s2  = BezierSegment_get__s2(seg);
    double cs1 = BezierSegment_get__cs1(seg);
    double cs2 = BezierSegment_get__cs2(seg);
    return 3*s2*t*t - 3*cs2*t*t + 6*cs2*(1 - t)*t - 6*cs1*(1 - t)*t - 3*s1*(1 - t)*(1 - t) + 3*cs1*(1 - t)*(1 - t);
}

/*gpufun*/
double BezierSegment_deriv_x(BezierSegment seg, double t){
    double x1  = BezierSegment_get__x1(seg);
    double x2  = BezierSegment_get__x2(seg);
    double cx1 = BezierSegment_get__cx1(seg);
    double cx2 = BezierSegment_get__cx2(seg);
    return 3*x2*t*t - 3*cx2*t*t + 6*cx2*(1 - t)*t - 6*cx1*(1 - t)*t - 3*x1*(1 - t)*(1 - t) + 3*cx1*(1 - t)*(1 - t);
}

/*gpufun*/
void BezierSegment_update_box(BezierSegment seg, BoundingBox box, double t1, double t2){
    double ts1 = BezierSegment_get__ts1(seg);
    double ts2 = BezierSegment_get__ts2(seg);
    double s1  = BezierSegment_func_s(seg, t1);
    double s2  = BezierSegment_func_s(seg, t2);
    double tx1 = BezierSegment_get__tx1(seg);
    double tx2 = BezierSegment_get__tx2(seg);
    double x1  = BezierSegment_func_x(seg, t1);
    double x2  = BezierSegment_func_x(seg, t2);
    double smin, smax, xmin, xmax;
    if (t1 <= ts1 && ts1 <= t2){
        double es1 = BezierSegment_get__es1(seg);
        if (t1 <= ts2 && ts2 <= t2){
            double es2 = BezierSegment_get__es2(seg);
            double s[4] = {s1, s2, es1, es2};
            sort_array_of_4_double(s);
            smin = s[0];
            smax = s[3];
        } else {
            double s[3] = {s1, s2, es1};
            sort_array_of_3_double(s);
            smin = s[0];
            smax = s[2];
        }
    } else if (t1 <= ts2 && ts2 <= t2){
        double es2 = BezierSegment_get__es2(seg);
        double s[3] = {s1, s2, es2};
        sort_array_of_3_double(s);
        smin = s[0];
        smax = s[2];
    } else {
        smin = MIN(s1, s2);
        smax = MAX(s1, s2);
    }
    if (t1 <= tx1 && tx1 <= t2){
        double ex1 = BezierSegment_get__ex1(seg);
        if (t1 <= tx2 && tx2 <= t2){
            double ex2 = BezierSegment_get__ex2(seg);
            double s[4] = {x1, x2, ex1, ex2};
            sort_array_of_4_double(s);
            xmin = s[0];
            xmax = s[3];
        } else {
            double s[3] = {x1, x2, ex1};
            sort_array_of_3_double(s);
            xmin = s[0];
            xmax = s[2];
        }
    } else if (t1 <= tx2 && tx2 <= t2){
        double ex2 = BezierSegment_get__ex2(seg);
        double s[3] = {x1, x2, ex2};
        sort_array_of_3_double(s);
        xmin = s[0];
        xmax = s[2];
    } else {
        xmin = MIN(x1, x2);
        xmax = MAX(x1, x2);
    }
    double rC = sqrt(smin*smin + smax*smax); // length of position vector to first vertex
    double sin_tC = xmin / rC; // angle of position vector to first vertex
    double cos_tC = smin / rC;
    double l = smax - smin;  // length of the box
    double w = xmax - xmin;  // width of the box
    double sin_tb = 0;       // orientation of the box (angle of length wrt horizontal)
    double cos_tb = 1;
    BoundingBox_set(box, rC, sin_tC, cos_tC, l, w, sin_tb, cos_tb);
}

/*gpufun*/
void BezierSegment_calculate_extrema(BezierSegment seg){
    double s1  = BezierSegment_get__s1(seg);
    double s2  = BezierSegment_get__s2(seg);
    double cs1 = BezierSegment_get__cs1(seg);
    double cs2 = BezierSegment_get__cs2(seg);
    double x1  = BezierSegment_get__x1(seg);
    double x2  = BezierSegment_get__x2(seg);
    double cx1 = BezierSegment_get__cx1(seg);
    double cx2 = BezierSegment_get__cx2(seg);
    double a = 3*(s2 - s1) - 9*(cs2 - cs1);
    double b = 6*(s1 + cs2 - 2*cs1);
    double c = -3*(s1 - cs1);
    double t1 = -1;
    double t2 = -1;
    if (fabs(a) < XC_GEOM_EPSILON){
        if (fabs(b) >= XC_GEOM_EPSILON){
            // Linear equation
            t1 = -c/b;
        } // else no extrema
    } else {
        // Quadratic equation
        double disc = b*b - 4*a*c;
        if (fabs(disc) < XC_GEOM_EPSILON){
            // One extremum
            t1 = -b/(2*a);
        } else if (disc >= XC_GEOM_EPSILON){
            // Two extrema
            t1 = (-b + sqrt(disc))/(2*a);
            t2 = (-b - sqrt(disc))/(2*a);
        } // else no extrema
    }
    BezierSegment_set__ts1(seg, t1);
    BezierSegment_set__ts2(seg, t2);
    BezierSegment_set__es1(seg, BezierSegment_func_s(seg, t1));
    BezierSegment_set__es2(seg, BezierSegment_func_s(seg, t2));
    a = 3*(x2 - x1) - 9*(cx2 - cx1);
    b = 6*(x1 + cx2 - 2*cx1);
    c = -3*(x1 - cx1);
    t1 = -1;
    t2 = -1;
    if (fabs(a) < XC_GEOM_EPSILON){
        if (fabs(b) >= XC_GEOM_EPSILON){
            // Linear equation
            t1 = -c/b;
        } // else no extrema
    } else {
        // Quadratic equation
        double disc = b*b - 4*a*c;
        if (fabs(disc) < XC_GEOM_EPSILON){
            // One extremum
            t1 = -b/(2*a);
        } else if (disc >= XC_GEOM_EPSILON){
            // Two extrema
            t1 = (-b + sqrt(disc))/(2*a);
            t2 = (-b - sqrt(disc))/(2*a);
        } // else no extrema
    }
    BezierSegment_set__tx1(seg, t1);
    BezierSegment_set__tx2(seg, t2);
    BezierSegment_set__ex1(seg, BezierSegment_func_x(seg, t1));
    BezierSegment_set__ex2(seg, BezierSegment_func_x(seg, t2));
}
// /*gpufun*/
// void _hit_s_bezier(BezierSegment seg, double t, double multiplicity, int8_t* n_hit, double* s){
//     double s1  = BezierSegment_get_s1(seg);
//     double s2  = BezierSegment_get_s2(seg);
//     double cs1 = BezierSegment_get_cs1(seg);
//     double cs2 = BezierSegment_get_cs2(seg);
//     double new_s = (1-t)*(1-t)*(1-t)*s1 + 3*(1-t)*(1-t)*t*cs1 + 3*(1-t)*t*t*cs2 + t*t*t*s2;
//     for (int8_t i = 0; i < multiplicity; i++) {
//         s[*n_hit] = new_s;
//         (*n_hit)++;
//     }
// }



#endif /* XCOLL_GEOM_TRAJ_BEZIER_H */