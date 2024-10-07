# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
from ..c_init import GeomCInit


# Line segment between two points (s1, x1) -- (s2, x2)
class LineSegment(xo.Struct):
    s1 = xo.Float64
    x1 = xo.Float64
    s2 = xo.Float64
    x2 = xo.Float64

    _depends_on = [GeomCInit]
    _extra_c_sources = [
        """
/*gpufun*/
void LineSegment_crossing_drift(LineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double m){
    // Get segment data
    double s1 = LineSegment_get_s1(seg);
    double x1 = LineSegment_get_x1(seg);
    double s2 = LineSegment_get_s2(seg);
    double x2 = LineSegment_get_x2(seg);
    double denom = (x2 - x1) - (s2 - s1)*m;
    if (fabs(denom) < XC_EPSILON){
        // Trajectory is parallel to the segment
        if (fabs((x0 - x1)/(s0 - s1) - m) < XC_EPSILON){
            // Trajectory overlaps with the segment
            // TODO: This is situational; we should return s1 if get_s_first and current_s if after current_s
            //       For now we hit twice (because we go nor IN nor OUT)
            s[*n_hit] = s1;
            (*n_hit)++;
            s[*n_hit] = s2;
            (*n_hit)++;
        } else {
            // No crossing
            return;
        }
    } else {
        double t = (x0 - x1 - (s0 - s1)*m) / denom;
        if (t >= 0 && t <= 1){
            s[*n_hit] = s1*(1-t) + s2*t;
            (*n_hit)++;
        }
    }
}
"""]

    def evaluate(self, t):
        s1 = self.s1
        x1 = self.x1
        s2 = self.s2
        x2 = self.x2
        t = np.array(t)
        mask = (t >= 0) & (t <= 1)
        return s1*(1-t[mask]) + s2*t[mask], x1*(1-t[mask]) + x2*t[mask]
