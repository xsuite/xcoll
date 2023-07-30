// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H
#include <math.h>
#include <stdio.h>



/*gpufun*/
void EverestCrystal_set_material(EverestCrystalData el, LocalParticle* part0){
    CrystalMaterialData material = EverestCrystalData_getp__material(el);
    RandomRutherfordData rng = EverestCrystalData_getp_rutherford_rng(el);
    RandomRutherford_set_by_xcoll_material(rng, (GeneralMaterialData) material);
}

/*gpufun*/
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t active      = EverestCrystalData_get_active(el);
    active            *= EverestCrystalData_get__tracking(el);
    double const inactive_front = EverestCrystalData_get_inactive_front(el);
    double const active_length  = EverestCrystalData_get_active_length(el);
    double const inactive_back  = EverestCrystalData_get_inactive_back(el);

    CrystalMaterialData material = EverestCrystalData_getp__material(el);
    RandomRutherfordData rng = EverestCrystalData_getp_rutherford_rng(el);

    // Crystal properties
    double const co_x       = EverestCrystalData_get_ref_x(el);
    double const co_y       = EverestCrystalData_get_ref_y(el);
    // TODO: use xtrack C-code for rotation element
    // TODO: we are ignoring the angle of the right jaw
    // TODO: is a crystal always one-sided...?
    double const sin_zL     = EverestCrystalData_get_sin_zL(el);
    double const cos_zL     = EverestCrystalData_get_cos_zL(el);
    double const sin_zR     = EverestCrystalData_get_sin_zR(el);
    double const cos_zR     = EverestCrystalData_get_cos_zR(el);
    if (fabs(sin_zL-sin_zR) > 1.e-10 || fabs(cos_zL-cos_zR) > 1.e-10 ){
        kill_all_particles(part0, XC_ERR_NOT_IMPLEMENTED);
    };
    double const c_aperture = EverestCrystalData_get_jaw_L(el) - EverestCrystalData_get_jaw_R(el);
    double const c_offset   = ( EverestCrystalData_get_jaw_L(el) + EverestCrystalData_get_jaw_R(el) ) /2;
    double const c_tilt0    = asin(EverestCrystalData_get_sin_yL(el));
    double const c_tilt1    = asin(EverestCrystalData_get_sin_yR(el));
    if (fabs(c_tilt1) > 1.e-10){
        kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
    };
    int    const side       = EverestCrystalData_get__side(el);
    double const bend_r     = EverestCrystalData_get__bending_radius(el);
    // TODO: cry_tilt should be given by jaw positions...?
    double const cry_tilt   = EverestCrystalData_get_align_angle(el) + c_tilt0;   // TODO: only left-sided crystals
    double const bend_ang   = EverestCrystalData_get__bending_angle(el);
    // double const cry_bend   = length/cry_rcurv; //final value (with corrected length)
    // THIS IS WRONG! Was a mistranslation from SixTrack 4 to SixTrack 5
    // Difference is so small that this was never caught.
    // Furthermore, we removed the adaptation of the scatter length, because
    // 1) it was implemented wrong (passed unnoticed due to small effect)
    // 2) we should not use the adapted scatter length, as we rotate the S-X frame, so
    //    we anyway have to drift the full length!    
    double const cry_alayer = EverestCrystalData_get_thick(el);
    double const xdim       = EverestCrystalData_get_xdim(el);
    double const ydim       = EverestCrystalData_get_ydim(el);
    double const cry_orient = EverestCrystalData_get__orient(el);
    double const cry_miscut = EverestCrystalData_get_miscut(el);

    // Impact table
    CollimatorImpactsData record = EverestCrystalData_getp_internal_record(el, part0);
    RecordIndex record_index = NULL;
    if (record){
        record_index = CollimatorImpactsData_getp__index(record);
    }

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, inactive_front + active_length + inactive_back);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
            int8_t rng_is_set  = assert_rng_set(part, RNG_ERR_SEEDS_NOT_SET);
            int8_t ruth_is_set = assert_rutherford_set(rng, part, RNG_ERR_RUTH_NOT_SET);

            if (is_tracking && rng_is_set && ruth_is_set) {
                // Drift inactive front
                Drift_single_particle(part, inactive_front);

                // Scatter

                // Move to collimator frame
                XYShift_single_particle(part, co_x, co_y);
                SRotation_single_particle(part, sin_zL, cos_zL);

                scatter_cry(part, active_length, material, rng, c_aperture, c_offset,
                            side, cry_tilt, bend_r, bend_ang, cry_alayer, xdim, ydim, cry_orient, 
                            cry_miscut, record, record_index);

                // Return from collimator frame
                SRotation_single_particle(part, -sin_zL, cos_zL);
                XYShift_single_particle(part, -co_x, -co_y);

                // Drift inactive back (only surviving particles)
                if (LocalParticle_get_state(part) > 0){
                    Drift_single_particle(part, inactive_back);
                }
            }
        }
    //end_per_particle_block
}


#endif /* XCOLL_EVEREST_CRYSTAL_H */
