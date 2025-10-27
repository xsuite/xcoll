// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_COLL_H
#define XCOLL_EVEREST_COLL_H
#include <math.h>
#include <stdio.h>


/*gpufun*/
int8_t EverestCollimatorData_get_record_impacts(EverestCollimatorData el){
    return EverestCollimatorData_get__record_interactions(el) % 2;
}

/*gpufun*/
int8_t EverestCollimatorData_get_record_exits(EverestCollimatorData el){
    return (EverestCollimatorData_get__record_interactions(el) >> 1) % 2;
}

/*gpufun*/
int8_t EverestCollimatorData_get_record_scatterings(EverestCollimatorData el){
    return (EverestCollimatorData_get__record_interactions(el) >> 2) % 2;
}

void EverestCollimator_set_material(EverestCollimatorData el){
    MaterialData material = EverestCollimatorData_getp__material(el);
    RandomRutherfordData rng = EverestCollimatorData_getp_rutherford_rng(el);
    RandomRutherford_set_by_xcoll_material(rng, material);
}


/*gpufun*/
CollimatorGeometry EverestCollimator_init_geometry(EverestCollimatorData el, LocalParticle* part0){
    CollimatorGeometry cg = (CollimatorGeometry) malloc(sizeof(CollimatorGeometry_));
    // Jaw corners (with tilts)
    cg->jaw_LU = EverestCollimatorData_get__jaw_LU(el);
    cg->jaw_RU = EverestCollimatorData_get__jaw_RU(el);
    // Get angles of jaws
    cg->sin_zL = EverestCollimatorData_get__sin_zL(el);
    cg->cos_zL = EverestCollimatorData_get__cos_zL(el);
    cg->sin_zR = EverestCollimatorData_get__sin_zR(el);
    cg->cos_zR = EverestCollimatorData_get__cos_zR(el);
    cg->sin_zDiff = EverestCollimatorData_get__sin_zDiff(el);
    cg->cos_zDiff = EverestCollimatorData_get__cos_zDiff(el);
    cg->jaws_parallel = EverestCollimatorData_get__jaws_parallel(el);
    // Tilts
    cg->sin_yL = EverestCollimatorData_get__sin_yL(el);
    cg->cos_yL = EverestCollimatorData_get__cos_yL(el);
    cg->sin_yR = EverestCollimatorData_get__sin_yR(el);
    cg->cos_yR = EverestCollimatorData_get__cos_yR(el);
    // Length and segments
    cg->length = EverestCollimatorData_get_length(el);
    cg->side   = EverestCollimatorData_get__side(el);
    double s_U, s_D, x_D;
    if (cg->side != -1){
        s_U = cg->length/2 * (1-cg->cos_yL);
        s_D = cg->length/2 * (1+cg->cos_yL);
        x_D = EverestCollimatorData_get__jaw_LD(el);
        cg->segments_L = create_jaw(s_U, cg->jaw_LU, s_D, x_D, cg->sin_yL/cg->cos_yL, 1);
    }
    if (cg->side != 1){
        s_U = cg->length/2 * (1-cg->cos_yR);
        s_D = cg->length/2 * (1+cg->cos_yR);
        x_D = EverestCollimatorData_get__jaw_RD(el);
        cg->segments_R = create_jaw(s_U, cg->jaw_RU, s_D, x_D, cg->sin_yR/cg->cos_yR, -1);
    }
    // Impact table
    cg->record = EverestCollimatorData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record){
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = EverestCollimatorData_get_record_impacts(el);
        cg->record_exits = EverestCollimatorData_get_record_exits(el);
    }
    return cg;
}

/*gpufun*/
void EverestCollimator_free(CollimatorGeometry restrict cg){
    if (cg->side != -1){
        destroy_jaw(cg->segments_L);
    }
    if (cg->side != 1){
        destroy_jaw(cg->segments_R);
    }
    free(cg);
}


// TODO: it would be great if we could set EverestData as an xofield, because then we could
// run this function at creation of the collimator instead of every turn
// Hmmmm this should be called whenever we change an xofield
/*gpufun*/
EverestCollData EverestCollimator_init(EverestCollimatorData el, LocalParticle* part0){
    EverestCollData coll = (EverestCollData) malloc(sizeof(EverestCollData_));
    // Random generator
    coll->rng = EverestCollimatorData_getp_rutherford_rng(el);
    // Impact table:  need it here to record interactions
    coll->record = EverestCollimatorData_getp_internal_record(el, part0);
    coll->record_index = NULL;
    coll->record_scatterings = 0;
    if (coll->record){
        coll->record_index = InteractionRecordData_getp__index(coll->record);
        coll->record_scatterings = EverestCollimatorData_get_record_scatterings(el);
    }
    coll->orient = 0;
    return coll;
}


/*gpufun*/
EverestData EverestCollimator_init_data(LocalParticle* part, MaterialData restrict material, EverestCollData coll){
    EverestData everest = (EverestData) malloc(sizeof(EverestData_));
    everest->coll = coll;
    everest->rescale_scattering = 1;
    // Preinitialise scattering parameters
    double energy = LocalParticle_get_energy(part) / 1e9; // energy in GeV
    calculate_scattering(everest, material, energy);
    calculate_ionisation_properties(everest, material, energy);
    return everest;
}


/*gpufun*/
void EverestCollimator_track_local_particle(EverestCollimatorData el, LocalParticle* part0) {
    int8_t active = EverestCollimatorData_get_active(el);
    active       *= EverestCollimatorData_get__tracking(el);
    double const length = EverestCollimatorData_get_length(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    EverestCollData coll;
    CollimatorGeometry cg;
    if (active){
        coll = EverestCollimator_init(el, part0);
        cg = EverestCollimator_init_geometry(el, part0);
    }

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                // Store s-location of start of collimator
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
                int8_t is_hit = hit_jaws_check_and_transform(part, cg);

                if (is_hit != 0) {
                    // Hit one of the jaws, so scatter
                    double remaining_length = length - LocalParticle_get_s(part);
                    // Scatter
                    MaterialData material = EverestCollimatorData_getp__material(el);
                    EverestData everest = EverestCollimator_init_data(part, material, coll);
                    pc_out = jaw(everest, material, part, pc_in, remaining_length, 1);
                    free(everest);
                }

                // Transform back to the lab frame
                hit_jaws_transform_back(is_hit, part, cg);
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
        EverestCollimator_free(cg);
        free(coll);
    }
}


#endif /* XCOLL_EVEREST_COLL_H */