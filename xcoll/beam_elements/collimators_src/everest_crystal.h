#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H
#include <math.h>


/*gpufun*/
void track_crystal(EverestCrystalData el, LocalParticle* part0) {

    // Material properties
    CrystalMaterialData material = EverestCrystalData_getp_material(el);
    double const zatom    = CrystalMaterialData_get_Z(material);
    double const emr      = CrystalMaterialData_get_nuclear_radius(material);
    double const hcut     = CrystalMaterialData_get_hcut(material);

    RandomGeneratorData rng = EverestCrystalData_getp_random_generator(el);
    RandomGeneratorData_set_rutherford(rng, zatom, emr, hcut);

    //start_per_particle_block (part0->part)
        int8_t is_tracking = xcoll_assert_tracking(part);
        int8_t rng_set     = xcoll_assert_rng_set(part);
        int8_t ruth_set    = xcoll_assert_rutherford_set(rng, part);

        if ( is_tracking && rng_set && ruth_set) {
            scatter_cry(el, part);
        }
    //end_per_particle_block
}

/*gpufun*/
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t is_active      = EverestCrystalData_get__active(el);
    is_active            *= EverestCrystalData_get__tracking(el);
    double const inactive_front = EverestCrystalData_get_inactive_front(el);
    double const active_length  = EverestCrystalData_get_active_length(el);
    double const inactive_back  = EverestCrystalData_get_inactive_back(el);

    if (!is_active){
        // Drift full length
        xcoll_drift_6d(part0, inactive_front+active_length+inactive_back);
    } else {
        xcoll_drift_6d(part0, inactive_front);
        track_crystal(el, part0);
        xcoll_drift_6d(part0, inactive_back);
    }
}

#endif
