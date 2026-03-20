// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_GLAUBER_SPLINES_H
#define XCOLL_EVEREST_GLAUBER_SPLINES_H
#include <math.h>


double Material_evaluate_glauber_spline(Material material, double sqrt_s) {
    double* knots = Material_get__cs_spline_knots(material);
    double* a     = Material_get__cs_spline_a(material);
    double* b     = Material_get__cs_spline_b(material);
    double* c     = Material_get__cs_spline_c(material);
    double* d     = Material_get__cs_spline_d(material);
    int n_points   = Material_get__n_points(material);
    double log_sqrt_s_min = Material_get__cs_log_sqrt_s_min(material);
    double log_step = Material_get__cs_log_step(material);

    // O(1) index
    int i = (int)((log(sqrt_s) - log_sqrt_s_min) / log_step);
    if (i < 0){
        i = 0;
    } else if (i >= n_points-1) {
        i = n_points - 2;
    }
    // Horner's method
    double dx = log(sqrt_s) - knots[i];
    return ((a[i]*dx + b[i])*dx + c[i])*dx + d[i];
}
s