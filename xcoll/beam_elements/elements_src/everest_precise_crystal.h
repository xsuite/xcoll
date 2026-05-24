// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_PRECISE_CRYSTAL_H
#define XCOLL_EVEREST_PRECISE_CRYSTAL_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU


/*gpufun*/
int8_t EverestPreciseCrystalData_get_record_impacts(EverestPreciseCrystalData el){
    return EverestPreciseCrystalData_get__record_interactions(el) % 2;
}

/*gpufun*/
int8_t EverestPreciseCrystalData_get_record_exits(EverestPreciseCrystalData el){
    return (EverestPreciseCrystalData_get__record_interactions(el) >> 1) % 2;
}

/*gpufun*/
int8_t EverestPreciseCrystalData_get_record_scatterings(EverestPreciseCrystalData el){
    return (EverestPreciseCrystalData_get__record_interactions(el) >> 2) % 2;
}


void EverestPreciseCrystal_set_material(EverestPreciseCrystalData el){
    MaterialData material = EverestPreciseCrystalData_getp__material(el);
    RandomRutherfordData rng = EverestPreciseCrystalData_getp_rutherford_rng(el);
    RandomRutherford_set_by_xcoll_material(rng, (MaterialData) material);
}


/*gpufun*/
CrystalGeometry EverestPreciseCrystal_init_geometry(EverestPreciseCrystalData el, LocalParticle* part0){
    CrystalGeometry cg = (CrystalGeometry) malloc(sizeof(CrystalGeometry_));
    cg->length = EverestPreciseCrystalData_get_length(el);
    cg->side   = EverestPreciseCrystalData_get__side(el);
    if (cg->side == 0){
        kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
        return cg;
    }
    double R   = EverestPreciseCrystalData_get__bending_radius(el);
    double t_R = EverestPreciseCrystalData_get__bending_angle(el);
    cg->bending_radius = R;
    cg->bending_angle  = t_R;
    cg->miscut_angle   = EverestPreciseCrystalData_get_miscut(el);
    cg->width  = EverestPreciseCrystalData_get__width(el);
    cg->height = EverestPreciseCrystalData_get__height(el);
    cg->jaw_U  = EverestPreciseCrystalData_get__jaw_U(el);
    cg->sin_z  = EverestPreciseCrystalData_get__sin_z(el);
    cg->cos_z  = EverestPreciseCrystalData_get__cos_z(el);
    cg->sin_y  = EverestPreciseCrystalData_get__sin_y(el);
    cg->cos_y  = EverestPreciseCrystalData_get__cos_y(el);
    // Segments
    if (cg->side == 1){
        cg->segments = create_crystal(cg->bending_radius, cg->width, cg->length, cg->jaw_U, \
                                        cg->sin_y, cg->cos_y);
    } else if (cg->side == -1){
        // jaw_U is the inner corner (shifted if right-sided crystal)
        cg->segments = create_crystal(cg->bending_radius, cg->width, cg->length, cg->jaw_U - cg->width, \
                                        cg->sin_y, cg->cos_y);
    }
    // // Jaw frame is always left-sided
    // cg->segments_jf = create_crystal(cg->bending_radius, cg->width, cg->length, 0, 0, 1);
    // Bend centre
    cg->s_B = 0;
    cg->x_B = cg->bending_radius;
    // Miscut centre
    cg->s_P = -R*sin(cg->miscut_angle);
    cg->x_P = R*cos(cg->miscut_angle);
    if (cg->side == 1 && R < 0){
        // If R<0, a left-sided crystal bends towards the beam
        cg->x_P = cg->x_P + cg->width;
        cg->x_B = cg->x_B + cg->width;
    } else if (cg->side == -1 && R > 0){
        // If R>0, a right-sided crystal bends towards the beam
        cg->x_P = cg->x_P - cg->width;
        cg->x_B = cg->x_B - cg->width;
    }
    if (cg->side == -1){
        // Mirror the crystal geometry
        cg->bending_radius = -cg->bending_radius;
        cg->bending_angle  = -cg->bending_angle;
        cg->miscut_angle   = -cg->miscut_angle;
        cg->x_P            = -cg->x_P;
        cg->x_B            = -cg->x_B;
    }
    // From here on, crystal geometry parameters can always be treated as left-sided.
    // Note that the segments are not mirrored, which is fine as get_s_of_first_crossing_with_vlimit
    // is absolute (not in the jaw reference frame). It is only after a hit is registered, that we
    // need to transform the particle to the jaw reference frame.
    double Rb;
    if (cg->miscut_angle > 0){
        Rb = R - cg->width;
    } else {
        Rb = R;
    }
    cg->t_VImax = atan( (Rb*sin(t_R) - cg->s_P) / (R - Rb*cos(t_R) - cg->x_P) );
    // Impact table
    cg->record = EverestPreciseCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record){
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = EverestPreciseCrystalData_get_record_impacts(el);
        cg->record_exits = EverestPreciseCrystalData_get_record_exits(el);
    }
    return cg;
}

/*gpufun*/
void EverestPreciseCrystal_free(CrystalGeometry restrict cg){
    destroy_crystal(cg->segments);
    free(cg);
}


// TODO: it would be great if we could set EverestData as an xofield, because then we could
// run this function at creation of the collimator instead of every turn
/*gpufun*/
EverestCollData EverestPreciseCrystal_init(
        EverestPreciseCrystalData el,
        MaterialData restrict material,
        LocalParticle* part0
){
    EverestCollData coll = (EverestCollData) malloc(sizeof(EverestCollData_));
    // Random generator
    coll->rng = EverestPreciseCrystalData_getp_rutherford_rng(el);
    // Impact table
    coll->record = EverestPreciseCrystalData_getp_internal_record(el, part0);
    coll->record_index = NULL;
    if (coll->record){
        coll->record_index = InteractionRecordData_getp__index(coll->record);
        coll->record_scatterings = EverestPreciseCrystalData_get_record_scatterings(el);
    }
    coll->orient = EverestPreciseCrystalData_get__orient(el);
    coll->is_precise = 1;
    // Configuration
    coll->method  = EverestPreciseCrystalData_get_method(el);
    coll->order   = EverestPreciseCrystalData_get_order(el);
    coll->variant = EverestPreciseCrystalData_get_variant(el);
    coll->n_steps = EverestPreciseCrystalData_get_n_steps(el);
    if (coll->n_steps <= 0){
        double _n_steps_auto = EverestPreciseCrystalData_get__n_steps_auto(el);
        double bpc = LocalParticle_get_p0c(part0)*LocalParticle_get_beta0(part0);
        coll->n_steps = (int64_t) (_n_steps_auto/sqrt(bpc));
   }
    // Parameters
    double aTF     = EverestPreciseCrystalData_get_aTF(el);
    double Umax    = EverestPreciseCrystalData_get_Umax(el);
    double alpha_i = EverestPreciseCrystalData_get_alpha_i(el);
    double beta_i  = EverestPreciseCrystalData_get_beta_i(el);
    double dp      = MaterialData_get__crystal_plane_distance(material)*0.002;  // [m]
    double eta     = MaterialData_get__eta(material);
    coll->x_c = sqrt(eta)*dp/2;
    coll->U_N = U_N_(Umax, dp, alpha_i, beta_i, beta_i/aTF);
    coll->aTF_over_beta = aTF / beta_i;
    coll->beta_over_aTF = beta_i / aTF;
    return coll;
}


/*gpufun*/
EverestData EverestPreciseCrystal_init_data(EverestCollData restrict coll){
    EverestData everest = (EverestData) malloc(sizeof(EverestData_));
    everest->coll = coll;
    everest->rescale_scattering = 1;
    everest->overbent = 0;
    return everest;
}


/*gpufun*/
void EverestPreciseCrystal_track_local_particle(EverestPreciseCrystalData el, LocalParticle* part0) {
    int8_t active = EverestPreciseCrystalData_get_active(el);
    active       *= EverestPreciseCrystalData_get__tracking(el);
    double length = EverestPreciseCrystalData_get_length(el);

    // Initialise collimator data
    EverestCollData coll;
    CrystalGeometry cg;
    MaterialData material;
    if (active){
        // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
        material = EverestPreciseCrystalData_getp__material(el);
        coll = EverestPreciseCrystal_init(el, material, part0);
        cg   = EverestPreciseCrystal_init_geometry(el, part0);

        // For info
        double const e0 = LocalParticle_get_energy0(part0);
        double t_c0  = _critical_angle0(material, e0);
        double Rcrit = _critical_radius(material, e0);
        double t_c = _critical_angle(coll, t_c0, Rcrit / fabs(cg->bending_radius));
        EverestPreciseCrystalData_set__critical_radius(el, Rcrit);
        EverestPreciseCrystalData_set__critical_angle(el, t_c);

        if (cg->width==0 || cg->height==0 || cg->bending_radius==0){
            kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
        }
    }

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                // Store s-location of start of crystal
                double const s_coll = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Store initial coordinates for updating later
                double const rvv_in  = LocalParticle_get_rvv(part);
#ifdef XCOLL_USE_EXACT
                double const xp_in   = LocalParticle_get_exact_xp(part);
                double const yp_in   = LocalParticle_get_exact_yp(part);
#else
                double const xp_in   = LocalParticle_get_xp(part);
                double const yp_in   = LocalParticle_get_yp(part);
#endif
                double const zeta_in = LocalParticle_get_zeta(part);
                double const p0c     = LocalParticle_get_p0c(part);
                double const delta   = LocalParticle_get_delta(part);
                double const qq0     = LocalParticle_get_charge_ratio(part);
                double const chi     = LocalParticle_get_chi(part);
                double const pc_in   = (1 + delta)*p0c*qq0/chi;
                double pc_out;

                // Check if hit on jaws
                int8_t is_hit = hit_crystal_check_and_transform(part, cg);

                if (is_hit != 0) {
                    // Hit one of the jaws, so scatter
                    double remaining_length = length - LocalParticle_get_s(part);
                    // Scatter
                    EverestData everest = EverestPreciseCrystal_init_data(coll);
                    pc_out = do_crystal_precise(everest, material, part, cg, pc_in/1.e9, remaining_length)*1.e9;
                    free(everest);
                }

                // Transform back to the lab frame
                hit_crystal_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);

                LocalParticle_set_zeta(part, zeta_in);

                // Hit and survived particles need correcting:
                if (is_hit!=0 && LocalParticle_get_state(part)>0){
                    double const rpp_old  = LocalParticle_get_rpp(part);
                    LocalParticle_update_delta(part, pc_out*chi/p0c/qq0 - 1);
                    // Keep angles constant (this is also correct for exact angles): px_new = px_old*(1 + δ_new)/(1 + δ_old)
                    double const scale = rpp_old / LocalParticle_get_rpp(part);
                    LocalParticle_scale_px(part, scale);
                    LocalParticle_scale_py(part, scale);

                    // Update zeta
#ifdef XCOLL_USE_EXACT
                    double xp  = LocalParticle_get_exact_xp(part);
                    double yp  = LocalParticle_get_exact_yp(part);
#else
                    double xp  = LocalParticle_get_xp(part);
                    double yp  = LocalParticle_get_yp(part);
#endif
                    double rvv = LocalParticle_get_rvv(part);
                    // First we drift half the length with the old angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, xp_in, yp_in, length/2) );
                    // then half the length with the new angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, xp, yp, length/2) );
                }
            }
        }
    //end_per_particle_block
    if (active){
        EverestPreciseCrystal_free(cg);
        free(coll);
    }
}


#endif /* XCOLL_EVEREST_PRECISE_CRYSTAL_H */
