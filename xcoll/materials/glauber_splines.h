// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_GLAUBER_SPLINES_H
#define XCOLL_EVEREST_GLAUBER_SPLINES_H
#include <math.h>

/*gpufun*/
double MaterialData_evaluate_glauber_spline(MaterialData material, double sqrt_s, int key) {
    // Key: 0 = total, 1 = production, 2 = elastic nucleus, 3 = elastic nucleon, 4 = single diffractive
    int n_points = MaterialData_get__n_points(material);
    double log_sqrt_s_min = MaterialData_get__cs_log_sqrt_s_min(material);
    double log_step = MaterialData_get__cs_log_step(material);

    if (n_points < 2 || log_step <= 0.0 || sqrt_s <= 0.0) {
        return 0.0;
    }

    // O(1) index in log(sqrt_s) space
    int i = (int)((log(sqrt_s) - log_sqrt_s_min) / log_step);
    if (i < 0) {
        i = 0;
    } else if (i >= n_points - 1) {
        i = n_points - 2;
    }

    double ai, bi, ci, di;
    switch (key) {
        case 0:
            ai = MaterialData_get__cs_tot_hA_a(material, i);
            bi = MaterialData_get__cs_tot_hA_b(material, i);
            ci = MaterialData_get__cs_tot_hA_c(material, i);
            di = MaterialData_get__cs_tot_hA_d(material, i);
            break;
        case 1:
            ai = MaterialData_get__cs_prod_hA_a(material, i);
            bi = MaterialData_get__cs_prod_hA_b(material, i);
            ci = MaterialData_get__cs_prod_hA_c(material, i);
            di = MaterialData_get__cs_prod_hA_d(material, i);
            break;
        case 2:
            ai = MaterialData_get__cs_el_hA_a(material, i);
            bi = MaterialData_get__cs_el_hA_b(material, i);
            ci = MaterialData_get__cs_el_hA_c(material, i);
            di = MaterialData_get__cs_el_hA_d(material, i);
            break;
        case 3:
            ai = MaterialData_get__cs_el_nucleon_a(material, i);
            bi = MaterialData_get__cs_el_nucleon_b(material, i);
            ci = MaterialData_get__cs_el_nucleon_c(material, i);
            di = MaterialData_get__cs_el_nucleon_d(material, i);
            break;
        case 4:
            ai = MaterialData_get__cs_sd_hA_a(material, i);
            bi = MaterialData_get__cs_sd_hA_b(material, i);
            ci = MaterialData_get__cs_sd_hA_c(material, i);
            di = MaterialData_get__cs_sd_hA_d(material, i);
            break;
        default:
            return 0.0; // Invalid key
    }

    double knot_i = MaterialData_get__cs_knots(material, i);
    double dx = log(sqrt_s) - knot_i;

    // Horner's method
    return ((ai*dx + bi)*dx + ci)*dx + di;
}
#endif /*XCOLL_EVEREST_GLAUBER_SPLINES_H*/