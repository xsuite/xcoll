// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_PRECISE_CHANNELLING_H
#define XCOLL_EVEREST_PRECISE_CHANNELLING_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU


/*gpufun*/
void calculate_potentials(EverestData restrict everest, CrystalGeometry restrict cg,
                          double pc){
    if (everest->coll->is_precise == 0) {
        return;
    }
    double bpc = pc*1e9;  // should actually be beta*pc
    double R = fabs(cg->bending_radius);
    double x_c = everest->coll->x_c;
    double aTF_over_beta = everest->coll->aTF_over_beta;
    double beta_over_aTF = everest->coll->beta_over_aTF;
    double U_N = everest->coll->U_N;
    double x_min = aTF_over_beta * asinh(-bpc / R / U_N * aTF_over_beta);
    everest->U_min  = U_simplemoliere(x_min, U_N, beta_over_aTF) + bpc/R * x_min;
    everest->U_well = U_simplemoliere(-x_c, U_N, beta_over_aTF) - bpc/R * x_c;
    if (everest->U_min >= everest->U_well) {
        everest->overbent = 1;
    }
}


/*gpufun*/
int8_t can_channel_precise(EverestData restrict everest, LocalParticle* part,
                           CrystalGeometry restrict cg, double x_local,
                           double xp_local, double pc) {
    if (everest->coll->is_precise == 0) {
        return 0;
    }
    if (everest->overbent) {
        return 0;
    }
    if (fabs(x_local) > everest->coll->x_c) {
        return 0;
    }
    double bpc = pc*1e9;  // should actually be beta*pc
    double R = fabs(cg->bending_radius);
    double U_N = everest->coll->U_N;
    double beta_over_aTF = everest->coll->beta_over_aTF;
    double U = U_simplemoliere(x_local, U_N, beta_over_aTF) + bpc/R*x_local;
    double E_T = E_T_simplemoliere(xp_local, bpc, U);
    if (E_T > everest->U_well || E_T < everest->U_min) {
        return 0;
    }
    return 1;
}


/*gpufun*/
double* calculate_local_coords(EverestData restrict everest, MaterialData restrict material,
                               LocalParticle* part, CrystalGeometry restrict cg) {
    double* result = (double*)malloc(3 * sizeof(double));

    double R = fabs(cg->bending_radius);
    double r = everest->r;
    double s = LocalParticle_get_s(part) - cg->s_P;
    double x = LocalParticle_get_x(part) - cg->x_P;

    // Precision is important! We are working with small differences of large numbers,
    // so we use FMA and HYPOT to ensure that we lose as little precision as possible.
    // The fma(x, y, z) function returns the result of x * y + z

    // First get r_local, the local x coordinate in the curved frame
    double delta = fma(-x, x, R * R);  // R^2 - x^2
    delta = fma(-s, s, delta);         // R^2 - s^2 - x^2 = (R - r)(R + r)
    double r_local = delta / (R + r);  // r_local = R - r

    // Then move to the channel frame, where x=0 is the middle of the channel
    // and positive x is towards the bending centre (higher potential)
    double dp = MaterialData_get__crystal_plane_distance(material)*0.002;  // [m]
    double channel_index = floor(r_local / dp);
    double x_local = fma(-channel_index - 0.5, dp, r_local);  // R - r - (channel_index + 0.5)*dp
#ifdef XCOLL_USE_EXACT
    double const xp = LocalParticle_get_exact_xp(part);
#else
    double const xp = LocalParticle_get_xp(part);
#endif
    double xp_local = xp - everest->t_I;
    double sign_R = copysign(1.0, R);
    xp_local = sign_R * xp_local;

    result[0] = x_local;
    result[1] = xp_local;
    result[2] = channel_index;
    return result;
}


/*gpufun*/
void set_lab_coords(MaterialData restrict material, LocalParticle* part,
                    CrystalGeometry restrict cg, double x_local, double xp_local,
                    double channel_index, double t_F) {
    double R = fabs(cg->bending_radius);
    double dp = MaterialData_get__crystal_plane_distance(material)*0.002;  // [m]
    double r_base = fma(-channel_index - 0.5, dp, R); // without small excursion x_local
    LocalParticle_set_x(part, r_base * cos(t_F) - x_local * cos(t_F) + cg->x_P);
    LocalParticle_set_s(part, r_base * sin(t_F) - x_local * sin(t_F) + cg->s_P);
    double sign_R = copysign(1.0, R);
#ifdef XCOLL_USE_EXACT
    LocalParticle_set_exact_xp(part, sign_R * xp_local + t_F);
#else
    LocalParticle_set_xp(part, sign_R * xp_local + t_F);
#endif
}


/*gpufun*/
void precise_channelling_transport(EverestData restrict everest, MaterialData restrict material,
                                   LocalParticle* part, CrystalGeometry restrict cg,
                                   double pc, double L_chan) {
    // Get parameters in the crystal frame
    double R = fabs(cg->bending_radius);
    double Ltot = everest->r * asin(cg->length/R);
    double bpc = pc*1e9;  // should actually be beta*pc
    double* coords = calculate_local_coords(everest, material, part, cg);
    double x = coords[0];
    double xp = coords[1];
    double channel_index = coords[2];
    free(coords);

    // Transport the particle to the end of channelling (to get y correct etc)
    // The distance from I to F is the chord length of the angle t_P: d = 2 r sin(t_P/2)
    // Hence the longitudinal distance (the length to be drifted) is the projection of this using the
    // xp at the start of channelling: s = 2 r sin(t_P/2)cos(t_P/2 + t_I)
    double t_I = everest->t_I;
    double t_P = L_chan / everest->r;
    double t_chord = t_I + t_P/2.;
    double drift_length = 2.*L_chan/t_P * sin(t_P/2.) * cos(t_chord);
    LocalParticle_set_xp(part, t_chord); // Angle at start of channelling
    Drift_single_particle_4d(part, drift_length);

    // Do exact transport in local frame
    int64_t n_steps = ceil(L_chan/Ltot * everest->coll->n_steps);
    const double ds = L_chan / (double)n_steps;
    double aTF_over_beta = everest->coll->aTF_over_beta;
    double beta_over_aTF = everest->coll->beta_over_aTF;
    double U_N = everest->coll->U_N;
    double U = U_simplemoliere(x, U_N, beta_over_aTF);

    for (int i = 0; i < n_steps; ++i) {
        if (everest->coll->method == 2) {
            fM2_apply_yoshida(
                ds, x, xp, bpc, R, aTF_over_beta, beta_over_aTF, U_N, U, sqrt(U),
                everest->coll->order, everest->coll->variant, &x, &xp
            );

        } else if (everest->coll->method == 3) {
            fM3_apply_yoshida(
                ds, x, xp, bpc, R, aTF_over_beta, beta_over_aTF, U_N, U, sqrt(U),
                everest->coll->order, everest->coll->variant, &x, &xp
            );

        } else if (everest->coll->method == 4) {
            fM4_apply_yoshida(
                ds, x, xp, bpc, R, beta_over_aTF, 2*aTF_over_beta,U_N,
                everest->coll->order, everest->coll->variant, &x, &xp
            );

        } else {
            kill_all_particles(part, XC_ERR_INVALID_XOFIELD);
            return;
        }
    }

    // And move the coordinates to the lab frame
    set_lab_coords(material, part, cg, x, xp, channel_index, t_I + t_P);
}

#endif /* XCOLL_EVEREST_PRECISE_CHANNELLING_H */
