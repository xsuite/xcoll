// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_H
#define XCOLL_EVEREST_H
#include <math.h>
#include <stdio.h>


/*gpufun*/
void EverestCollimator_track_local_particle(EverestCollimatorData el, LocalParticle* part0) {
    int8_t active = EverestCollimatorData_get_active(el);
    active       *= EverestCollimatorData_get__tracking(el);
    double const inactive_front = EverestCollimatorData_get_inactive_front(el);
    double const active_length  = EverestCollimatorData_get_active_length(el);
    double const inactive_back  = EverestCollimatorData_get_inactive_back(el);

    MaterialData material = EverestCollimatorData_getp__material(el);
    RandomRutherfordData rng = EverestCollimatorData_getp_rutherford_rng(el);
    RandomRutherford_set_by_xcoll_material(rng, (GeneralMaterialData) material);

    // Collimator properties
    double const length     = EverestCollimatorData_get_active_length(el);
    double const co_x       = EverestCollimatorData_get_ref_x(el);
    double const co_y       = EverestCollimatorData_get_ref_y(el);
    // TODO: use xtrack C-code for rotation element
    // TODO: we are ignoring the angle of the right jaw
    double const sin_zL     = EverestCollimatorData_get_sin_zL(el);
    double const cos_zL     = EverestCollimatorData_get_cos_zL(el);
    double const sin_zR     = EverestCollimatorData_get_sin_zR(el);
    double const cos_zR     = EverestCollimatorData_get_cos_zR(el);
    if (abs(sin_zL-sin_zR) > 1.e-10 || abs(cos_zL-cos_zR) > 1.e-10 ){
        kill_all_particles(part0, XC_ERR_NOT_IMPLEMENTED);
    };
    double const c_aperture = EverestCollimatorData_get_jaw_L(el) - EverestCollimatorData_get_jaw_R(el);
    double const c_offset   = ( EverestCollimatorData_get_jaw_L(el) + EverestCollimatorData_get_jaw_R(el) ) /2;
    double const c_tilt0    = asin(EverestCollimatorData_get_sin_yL(el));
    double const c_tilt1    = asin(EverestCollimatorData_get_sin_yR(el));
    int    const side       = EverestCollimatorData_get__side(el);

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, inactive_front+active_length+inactive_back);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
            int8_t rng_set     = assert_rng_set(part, RNG_ERR_SEEDS_NOT_SET);
            int8_t ruth_set    = assert_rutherford_set(rng, part, RNG_ERR_RUTH_NOT_SET);

            if (is_tracking && rng_set && ruth_set) {
                // Drift inactive front
                Drift_single_particle(part, inactive_front);

                // Scattering parameters
                double const energy0 = LocalParticle_get_energy0(part) / 1e9; // Reference energy in GeV
                struct ScatteringParameters scat = calculate_scattering(energy0, material);

                // Move to closed orbit
                LocalParticle_add_to_x(part, -co_x);
                LocalParticle_add_to_y(part, -co_y);
                
                scatter(part, length, material, rng, scat, cos_zL, sin_zL, c_aperture, c_offset, c_tilt0, c_tilt1, side);

                // Return from closed orbit
                LocalParticle_add_to_x(part, co_x);
                LocalParticle_add_to_y(part, co_y);

                // Drift inactive back (only surviving particles)
                if (LocalParticle_get_state(part) > 0){
                    Drift_single_particle(part, inactive_back);
                }
            }
        }
    //end_per_particle_block
}


#endif /* XCOLL_EVEREST_H */
