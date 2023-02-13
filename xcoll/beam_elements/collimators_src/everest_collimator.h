#ifndef XCOLL_EVEREST_H
#define XCOLL_EVEREST_H
#include <math.h>
#include <stdio.h>



/*gpufun*/
void track_collimator(EverestCollimatorData el, LocalParticle* part0) {

    // Material properties
    MaterialData material = EverestCollimatorData_getp_material(el);
    double const zatom    = MaterialData_get_Z(material);
    double const emr      = MaterialData_get_nuclear_radius(material);
    double const hcut     = MaterialData_get_hcut(material);

    RandomGeneratorData rng = EverestCollimatorData_getp_random_generator(el);
    RandomGeneratorData_set_rutherford(rng, zatom, emr, hcut);

    //start_per_particle_block (part0->part)
        int8_t is_tracking = xcoll_assert_tracking(part);
        int8_t rng_set     = xcoll_assert_rng_set(part);
        int8_t ruth_set    = xcoll_assert_rutherford_set(rng, part);

        if ( is_tracking && rng_set && ruth_set) {
            // Calculate scattering parameters
            double const energy0 = LocalParticle_get_energy0(part) / 1e9; // Reference energy in GeV
            struct ScatteringParameters scat = calculate_scattering(energy0, material);
            scatter(el, part, scat);
        }
    //end_per_particle_block
}

/*gpufun*/
void EverestCollimator_track_local_particle(EverestCollimatorData el, LocalParticle* part0) {
    int8_t is_active = EverestCollimatorData_get__active(el);
    is_active       *= EverestCollimatorData_get__tracking(el);

    double const inactive_front = EverestCollimatorData_get_inactive_front(el);
    double const active_length  = EverestCollimatorData_get_active_length(el);
    double const inactive_back  = EverestCollimatorData_get_inactive_back(el);

    if (!is_active){
        // Drift full length
        xcoll_drift_6d(part0, inactive_front+active_length+inactive_back);
    } else {
        xcoll_drift_6d(part0, inactive_front);
        track_collimator(el, part0);
        //TODO: check that this only drifts surviving particles!!!
        xcoll_drift_6d(part0, inactive_back);
    }
}

#endif
