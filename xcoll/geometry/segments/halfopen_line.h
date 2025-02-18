// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SEG_HALFOPENLINE_H
#define XCOLL_GEOM_SEG_HALFOPENLINE_H
#define XC_HALFOPENLINE_CROSSINGS 2


/*gpufun*/
double HalfOpenLineSegment_func_s(HalfOpenLineSegment seg, double t){
    double s1 = HalfOpenLineSegment_get_s1(seg);
    double cos_t1 = HalfOpenLineSegment_get_cos_t1(seg);
    return s1 + t*cos_t1;
}

/*gpufun*/
double HalfOpenLineSegment_func_x(HalfOpenLineSegment seg, double t){
    double x1 = HalfOpenLineSegment_get_x1(seg);
    double sin_t1 = HalfOpenLineSegment_get_sin_t1(seg);
    return x1 + t*sin_t1;
}

/*gpufun*/
double HalfOpenLineSegment_deriv_s(HalfOpenLineSegment seg, double t){
    UNUSED(t);
    return HalfOpenLineSegment_get_cos_t1(seg);
}

/*gpufun*/
double HalfOpenLineSegment_deriv_x(HalfOpenLineSegment seg, double t){
    UNUSED(t);
    return HalfOpenLineSegment_get_sin_t1(seg);
}

/*gpufun*/
int8_t HalfOpenLineSegment_func_in_domain(HalfOpenLineSegment seg, double t){
    return t >= 0;
}





/*gpufun*/
void HalfOpenLineSegment_crossing_drift(HalfOpenLineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
    // Get segment data
    double s1 = HalfOpenLineSegment_get_s(seg);
    double x1 = HalfOpenLineSegment_get_x(seg);
    double s2 = s1 + cos(HalfOpenLineSegment_get_t(seg));
    double x2 = x1 + sin(HalfOpenLineSegment_get_t(seg));
    double denom = (x2 - x1) - (s2 - s1)*xm;
    if (fabs(denom) < XC_EPSILON){
        // Trajectory is parallel to the segment
        if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_EPSILON){
            // Trajectory overlaps with the segment
            // TODO: This is situational; we should return s1 if get_s_first and current_s if after current_s
            //       For now we hit twice (because we go nor IN nor OUT)
            s[*n_hit] = s1;
            (*n_hit)++;
            s[*n_hit] = s1;
            (*n_hit)++;
        } else {
            // No hit
            return;
        }
    } else {
        double t = (x0 - x1 - (s0 - s1)*xm) / denom;
        if (t >= 0){  // We do not check for t<=1 as it is a half-open segment
            s[*n_hit] = s1*(1-t) + s2*t;
            (*n_hit)++;
        }
    }
}

/*gpufun*/
void HalfOpenLineSegment_crossing_mcs(HalfOpenLineSegment seg, int8_t* n_hit, double* s, const double* Ax, const double Xo){
   return grid_search_and_newton()
}  
#endif /* XCOLL_GEOM_SEG_HALFOPENLINE_H */