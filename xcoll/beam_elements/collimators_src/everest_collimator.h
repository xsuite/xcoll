#ifndef XCOLL_EVEREST_H
#define XCOLL_EVEREST_H
#include <math.h>
#include <stdio.h>



/*gpufun*/
void EverestCollimator_track_local_particle(EverestCollimatorData el, LocalParticle* part0) {
    int8_t is_active = EverestCollimatorData_get__active(el);
    is_active       *= EverestCollimatorData_get__tracking(el);
    double const inactive_front = EverestCollimatorData_get_inactive_front(el);
    double const active_length  = EverestCollimatorData_get_active_length(el);
    double const inactive_back  = EverestCollimatorData_get_inactive_back(el);

    MaterialData material = EverestCollimatorData_getp_material(el);
    RandomGeneratorData rng = EverestCollimatorData_getp_random_generator(el);
    RandomGeneratorData_set_rutherford_by_xcoll_material(rng, (GeneralMaterialData) material);

    //start_per_particle_block (part0->part)
        if (!is_active){
            // Drift full length
            Drift_single_particle(part, inactive_front+active_length+inactive_back);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = xcoll_assert_tracking(part);
            int8_t rng_set     = xcoll_assert_rng_set(part);
            int8_t ruth_set    = xcoll_assert_rutherford_set(rng, part);

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

#endif
