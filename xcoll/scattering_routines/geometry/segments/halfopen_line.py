# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
from ..c_init import GeomCInit


# Half-open line segment from a point (s, x) to +/-infinity along an angle
class HalfOpenLineSegment(xo.Struct):
    s = xo.Float64
    x = xo.Float64
    t = xo.Float64 # angle (wrt s-axis) towards inf

    _depends_on = [GeomCInit]
    _extra_c_sources = [
        """
/*gpufun*/
void HalfOpenLineSegment_crossing_drift(HalfOpenLineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double m){
    // Get segment data
    double s1 = HalfOpenLineSegment_get_s(seg);
    double x1 = HalfOpenLineSegment_get_x(seg);
    double s2 = s1 + cos(HalfOpenLineSegment_get_t(seg));
    double x2 = x1 + sin(HalfOpenLineSegment_get_t(seg));
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
            // No hit
            return;
        }
    } else {
        double t = (x0 - x1 - (s0 - s1)*m) / denom;
        if (t >= 0){  // We do not check for t<=1 as it is a half-open segment
            s[*n_hit] = s1*(1-t) + s2*t;
            (*n_hit)++;
        }
    }
}
"""]

    def evaluate(self, t):
        s1 = self.s
        x1 = self.x
        s2 = s1 + np.cos(self.t)
        x2 = x1 + np.sin(self.t)
        t = np.array(t)
        mask = t >= 0
        return s1*(1-t[mask]) + s2*t[mask], x1*(1-t[mask]) + x2*t[mask]
