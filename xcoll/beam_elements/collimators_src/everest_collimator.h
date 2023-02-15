// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_H
#define XCOLL_EVEREST_H
#include <math.h>
#include <stdio.h>

/*gpukern*/
void RandomRutherfordData_set_by_xcoll_material(RandomRutherfordData ran, GeneralMaterialData material){
    double const zatom    = GeneralMaterialData_get_Z(material);
    double const emr      = GeneralMaterialData_get_nuclear_radius(material);
    double const hcut     = GeneralMaterialData_get_hcut(material);
    double const lcut     = 0.0009982;
    double const c = 0.8561e3; // TODO: Where tha fuck does this come from??
    double A = pow(zatom,2);
    double B = c*pow(emr,2);
    RandomRutherfordData_set(ran, A, B, lcut, hcut);
}

/*gpufun*/
void EverestCollimator_track_local_particle(EverestCollimatorData el, LocalParticle* part0) {
    int8_t is_active = EverestCollimatorData_get__active(el);
    is_active       *= EverestCollimatorData_get__tracking(el);
    double const inactive_front = EverestCollimatorData_get_inactive_front(el);
    double const active_length  = EverestCollimatorData_get_active_length(el);
    double const inactive_back  = EverestCollimatorData_get_inactive_back(el);

    MaterialData material = EverestCollimatorData_getp_material(el);
    RandomRutherfordData rng = EverestCollimatorData_getp_rutherford_rng(el);
    RandomRutherfordData_set_by_xcoll_material(rng, (GeneralMaterialData) material);

    //start_per_particle_block (part0->part)
        if (!is_active){
            // Drift full length
            Drift_single_particle(part, inactive_front+active_length+inactive_back);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
            int8_t rng_set     = assert_rng_set(part, RNG_ERR_SEEDS_NOT_SET);
            int8_t ruth_set    = assert_rutherford_set(rng, part, RNG_ERR_RUTH_NOT_SET);

            if ( is_tracking && rng_set && ruth_set) {
                // Drift inactive front
                Drift_single_particle(part, inactive_front);

                // Scatter
                double const energy0 = LocalParticle_get_energy0(part) / 1e9; // Reference energy in GeV
                struct ScatteringParameters scat = calculate_scattering(energy0, material);
                scatter(el, part, scat);

                // Drift inactive back (only surviving particles)
                if (LocalParticle_get_state(part) > 0){
                    Drift_single_particle(part, inactive_back);
                }
            }
        }
    //end_per_particle_block
}


#endif /* XCOLL_EVEREST_H */
