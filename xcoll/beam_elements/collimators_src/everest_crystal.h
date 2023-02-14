// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H
#include <math.h>
#include <stdio.h>


/*gpufun*/
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t is_active      = EverestCrystalData_get__active(el);
    is_active            *= EverestCrystalData_get__tracking(el);
    double const inactive_front = EverestCrystalData_get_inactive_front(el);
    double const active_length  = EverestCrystalData_get_active_length(el);
    double const inactive_back  = EverestCrystalData_get_inactive_back(el);

    CrystalMaterialData material = EverestCrystalData_getp_material(el);
    RandomRutherfordData rng = EverestCrystalData_getp_rutherford_rng(el);
    RandomRutherfordData_set_by_xcoll_material(rng, (GeneralMaterialData) material);

    //start_per_particle_block (part0->part)
        if (!is_active){
            // Drift full length
            Drift_single_particle(part, inactive_front+active_length+inactive_back);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, xcoll_state_invalid_tracking);
            int8_t rng_set     = assert_rng_set(part, xcoll_state_rng_seeds_not_set);
            int8_t ruth_set    = assert_rutherford_set(rng, part, xcoll_state_rng_rutherford_not_set);

            if ( is_tracking && rng_set && ruth_set) {
                // Drift inactive front
                Drift_single_particle(part, inactive_front);

                // Scatter
                scatter_cry(el, part);

                // Drift inactive back (only surviving particles)
                if (LocalParticle_get_state(part) > 0){
                    Drift_single_particle(part, inactive_back);
                }
            }
        }
    //end_per_particle_block
}


#endif /* XCOLL_EVEREST_CRYSTAL_H */
