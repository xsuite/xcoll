// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_GLAUBER_SPLINES_H
#define XCOLL_EVEREST_GLAUBER_SPLINES_H
#include <math.h>

/*gpufun*/
double MaterialData_evaluate_glauber_spline(MaterialData material, double sqrt_s, int key, int particle_id) {
    // Key: 0 = total, 1 = production, 2 = elastic nucleus, 3 = elastic nucleon, 4 = single diffractive
    double ai, bi, ci, di, knot_i, log_sqrt_s_min, log_step, n_points;
    int i;
    switch (particle_id) {
        case 0: // proton
            n_points = MaterialData_get__n_points_pp(material);
            log_sqrt_s_min = MaterialData_get__cs_log_sqrt_s_min_pp(material);
            log_step = MaterialData_get__cs_log_step_pp(material);

            if (n_points < 2 || log_step <= 0.0 || sqrt_s <= 0.0) {
                return 0.0;
            }
            i = (int)((log(sqrt_s) - log_sqrt_s_min) / log_step);
            if (i < 0) {
                i = 0;
            } else if (i >= n_points - 1) {
                i = n_points - 2;
            }

            knot_i = MaterialData_get__cs_knots_pp(material, i);
            switch (key) {
                case 0:
                    ai = MaterialData_get__cs_tot_pp_a(material, i);
                    bi = MaterialData_get__cs_tot_pp_b(material, i);
                    ci = MaterialData_get__cs_tot_pp_c(material, i);
                    di = MaterialData_get__cs_tot_pp_d(material, i);
                    break;
                case 1:
                    ai = MaterialData_get__cs_prod_pp_a(material, i);
                    bi = MaterialData_get__cs_prod_pp_b(material, i);
                    ci = MaterialData_get__cs_prod_pp_c(material, i);
                    di = MaterialData_get__cs_prod_pp_d(material, i);
                    break;
                case 2:
                    ai = MaterialData_get__cs_el_pp_a(material, i);
                    bi = MaterialData_get__cs_el_pp_b(material, i);
                    ci = MaterialData_get__cs_el_pp_c(material, i);
                    di = MaterialData_get__cs_el_pp_d(material, i);
                    break;
                case 3:
                    ai = MaterialData_get__cs_el_nucleon_pp_a(material, i);
                    bi = MaterialData_get__cs_el_nucleon_pp_b(material, i);
                    ci = MaterialData_get__cs_el_nucleon_pp_c(material, i);
                    di = MaterialData_get__cs_el_nucleon_pp_d(material, i);
                    break;
                case 4:
                    ai = MaterialData_get__cs_sd_pp_a(material, i);
                    bi = MaterialData_get__cs_sd_pp_b(material, i);
                    ci = MaterialData_get__cs_sd_pp_c(material, i);
                    di = MaterialData_get__cs_sd_pp_d(material, i);
                    break;
                default:
                    return 0.0; // Invalid key
            }
            break;
        case 1: // Kaon minus
            n_points = MaterialData_get__n_points_kmin(material);
            log_sqrt_s_min = MaterialData_get__cs_log_sqrt_s_min_kmin(material);
            log_step = MaterialData_get__cs_log_step_kmin(material);

            if (n_points < 2 || log_step <= 0.0 || sqrt_s <= 0.0) {
                return 0.0;
            }
            i = (int)((log(sqrt_s) - log_sqrt_s_min) / log_step);
            if (i < 0) {
                i = 0;
            } else if (i >= n_points - 1) {
                i = n_points - 2;
            }
            knot_i = MaterialData_get__cs_knots_kmin(material, i);
            switch (key) {
                case 0:
                    ai = MaterialData_get__cs_tot_kmin_a(material, i);
                    bi = MaterialData_get__cs_tot_kmin_b(material, i);
                    ci = MaterialData_get__cs_tot_kmin_c(material, i);
                    di = MaterialData_get__cs_tot_kmin_d(material, i);
                    break;
                case 1:
                    ai = MaterialData_get__cs_prod_kmin_a(material, i);
                    bi = MaterialData_get__cs_prod_kmin_b(material, i);
                    ci = MaterialData_get__cs_prod_kmin_c(material, i);
                    di = MaterialData_get__cs_prod_kmin_d(material, i);
                    break;
                case 2:
                    ai = MaterialData_get__cs_el_kmin_a(material, i);
                    bi = MaterialData_get__cs_el_kmin_b(material, i);
                    ci = MaterialData_get__cs_el_kmin_c(material, i);
                    di = MaterialData_get__cs_el_kmin_d(material, i);
                    break;
                case 3:
                    ai = MaterialData_get__cs_el_nucleon_kmin_a(material, i);
                    bi = MaterialData_get__cs_el_nucleon_kmin_b(material, i);
                    ci = MaterialData_get__cs_el_nucleon_kmin_c(material, i);
                    di = MaterialData_get__cs_el_nucleon_kmin_d(material, i);
                    break;
                case 4:
                    ai = MaterialData_get__cs_sd_kmin_a(material, i);
                    bi = MaterialData_get__cs_sd_kmin_b(material, i);
                    ci = MaterialData_get__cs_sd_kmin_c(material, i);
                    di = MaterialData_get__cs_sd_kmin_d(material, i);
                    break;
                default:
                    return 0.0; // Invalid key
            }
            break;
        case 2: // Kaon plus
            n_points = MaterialData_get__n_points_kplus(material);
            log_sqrt_s_min = MaterialData_get__cs_log_sqrt_s_min_kplus(material);
            log_step = MaterialData_get__cs_log_step_kplus(material);

            if (n_points < 2 || log_step <= 0.0 || sqrt_s <= 0.0) {
                return 0.0;
            }
            i = (int)((log(sqrt_s) - log_sqrt_s_min) / log_step);
            if (i < 0) {
                i = 0;
            } else if (i >= n_points - 1) {
                i = n_points - 2;
            }
            knot_i = MaterialData_get__cs_knots_kplus(material, i);
            
            switch (key) {
                case 0:
                    ai = MaterialData_get__cs_tot_kplus_a(material, i);
                    bi = MaterialData_get__cs_tot_kplus_b(material, i);
                    ci = MaterialData_get__cs_tot_kplus_c(material, i);
                    di = MaterialData_get__cs_tot_kplus_d(material, i);
                    break;
                case 1:
                    ai = MaterialData_get__cs_prod_kplus_a(material, i);
                    bi = MaterialData_get__cs_prod_kplus_b(material, i);
                    ci = MaterialData_get__cs_prod_kplus_c(material, i);
                    di = MaterialData_get__cs_prod_kplus_d(material, i);
                    break;
                case 2:
                    ai = MaterialData_get__cs_el_kplus_a(material, i);
                    bi = MaterialData_get__cs_el_kplus_b(material, i);
                    ci = MaterialData_get__cs_el_kplus_c(material, i);
                    di = MaterialData_get__cs_el_kplus_d(material, i);
                    break;
                case 3:
                    ai = MaterialData_get__cs_el_nucleon_kplus_a(material, i);
                    bi = MaterialData_get__cs_el_nucleon_kplus_b(material, i);
                    ci = MaterialData_get__cs_el_nucleon_kplus_c(material, i);
                    di = MaterialData_get__cs_el_nucleon_kplus_d(material, i);
                    break;
                case 4:
                    ai = MaterialData_get__cs_sd_kplus_a(material, i);
                    bi = MaterialData_get__cs_sd_kplus_b(material, i);
                    ci = MaterialData_get__cs_sd_kplus_c(material, i);
                    di = MaterialData_get__cs_sd_kplus_d(material, i);
                    break;
                default:
                    return 0.0; // Invalid key
            }
            break;
        case 3: // Pion minus
            n_points = MaterialData_get__n_points_pimin(material);
            log_sqrt_s_min = MaterialData_get__cs_log_sqrt_s_min_pimin(material);
            log_step = MaterialData_get__cs_log_step_pimin(material);

            if (n_points < 2 || log_step <= 0.0 || sqrt_s <= 0.0) {
                return 0.0;
            }
            i = (int)((log(sqrt_s) - log_sqrt_s_min) / log_step);
            if (i < 0) {
                i = 0;
            } else if (i >= n_points - 1) {
                i = n_points - 2;
            } 
            knot_i = MaterialData_get__cs_knots_pimin(material, i);
            switch (key) {
                case 0:
                    ai = MaterialData_get__cs_tot_pimin_a(material, i);
                    bi = MaterialData_get__cs_tot_pimin_b(material, i);
                    ci = MaterialData_get__cs_tot_pimin_c(material, i);
                    di = MaterialData_get__cs_tot_pimin_d(material, i);
                    break;
                case 1:
                    ai = MaterialData_get__cs_prod_pimin_a(material, i);
                    bi = MaterialData_get__cs_prod_pimin_b(material, i);
                    ci = MaterialData_get__cs_prod_pimin_c(material, i);
                    di = MaterialData_get__cs_prod_pimin_d(material, i);
                    break;
                case 2:
                    ai = MaterialData_get__cs_el_pimin_a(material, i);
                    bi = MaterialData_get__cs_el_pimin_b(material, i);
                    ci = MaterialData_get__cs_el_pimin_c(material, i);
                    di = MaterialData_get__cs_el_pimin_d(material, i);
                    break;
                case 3:
                    ai = MaterialData_get__cs_el_nucleon_pimin_a(material, i);
                    bi = MaterialData_get__cs_el_nucleon_pimin_b(material, i);
                    ci = MaterialData_get__cs_el_nucleon_pimin_c(material, i);
                    di = MaterialData_get__cs_el_nucleon_pimin_d(material, i);
                    break;
                case 4:
                    ai = MaterialData_get__cs_sd_pimin_a(material, i);
                    bi = MaterialData_get__cs_sd_pimin_b(material, i);
                    ci = MaterialData_get__cs_sd_pimin_c(material, i);
                    di = MaterialData_get__cs_sd_pimin_d(material, i);
                    break;
                default:
                    return 0.0; // Invalid key
            }
             break;
        case 4: // Pion plus
            n_points = MaterialData_get__n_points_piplus(material);
            log_sqrt_s_min = MaterialData_get__cs_log_sqrt_s_min_piplus(material);
            log_step = MaterialData_get__cs_log_step_piplus(material);

            if (n_points < 2 || log_step <= 0.0 || sqrt_s <= 0.0) {
                return 0.0;
            }
            i = (int)((log(sqrt_s) - log_sqrt_s_min) / log_step);
            if (i < 0) {
                i = 0;
            } else if (i >= n_points - 1) {
                i = n_points - 2;
            }
            knot_i = MaterialData_get__cs_knots_piplus(material, i);
            switch (key) {
                case 0:
                    ai = MaterialData_get__cs_tot_piplus_a(material, i);
                    bi = MaterialData_get__cs_tot_piplus_b(material, i);
                    ci = MaterialData_get__cs_tot_piplus_c(material, i);
                    di = MaterialData_get__cs_tot_piplus_d(material, i);
                    break;
                case 1:
                    ai = MaterialData_get__cs_prod_piplus_a(material, i);
                    bi = MaterialData_get__cs_prod_piplus_b(material, i);
                    ci = MaterialData_get__cs_prod_piplus_c(material, i);
                    di = MaterialData_get__cs_prod_piplus_d(material, i);
                    break;
                case 2:
                    ai = MaterialData_get__cs_el_piplus_a(material, i);
                    bi = MaterialData_get__cs_el_piplus_b(material, i);
                    ci = MaterialData_get__cs_el_piplus_c(material, i);
                    di = MaterialData_get__cs_el_piplus_d(material, i);
                    break;
                case 3:
                    ai = MaterialData_get__cs_el_nucleon_piplus_a(material, i);
                    bi = MaterialData_get__cs_el_nucleon_piplus_b(material, i);
                    ci = MaterialData_get__cs_el_nucleon_piplus_c(material, i);
                    di = MaterialData_get__cs_el_nucleon_piplus_d(material, i);
                    break;
                case 4:
                    ai = MaterialData_get__cs_sd_piplus_a(material, i);
                    bi = MaterialData_get__cs_sd_piplus_b(material, i);
                    ci = MaterialData_get__cs_sd_piplus_c(material, i);
                    di = MaterialData_get__cs_sd_piplus_d(material, i);
                    break;
                default:
                    return 0.0; // Invalid key
            }
            break;
        default:
            return 0.0; // Invalid particle_id
    }

    double dx = log(sqrt_s) - knot_i;

    // Horner's method
    return ((ai*dx + bi)*dx + ci)*dx + di;
}
#endif /*XCOLL_EVEREST_GLAUBER_SPLINES_H*/