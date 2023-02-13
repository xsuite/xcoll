#ifndef XCOLL_EVEREST_H
#define XCOLL_EVEREST_H
#include <math.h>
#include <stdio.h>


// TODO:
//    Do not split 4d and zeta in drifts
//    Use drift function from xtrack Drift element (call its C function)
//    Use rotation function from xtrack XYRotation element (call its C function)

/*gpufun*/
void drift_6d(LocalParticle* part0, double length) {
    //start_per_particle_block (part0->part)
        double const rpp    = LocalParticle_get_rpp(part);
        double const rv0v   = 1./LocalParticle_get_rvv(part);
        double const xp     = LocalParticle_get_px(part) * rpp;
        double const yp     = LocalParticle_get_py(part) * rpp;
        double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

        LocalParticle_add_to_x(part, xp * length );
        LocalParticle_add_to_y(part, yp * length );
        LocalParticle_add_to_s(part, length);
        LocalParticle_add_to_zeta(part, length * dzeta );
    //end_per_particle_block
}

// TODO: 
// Write impacts



/*gpufun*/
void track_collimator(EverestCollimatorData el, LocalParticle* part0) {

    // Material properties
    MaterialData material = EverestCollimatorData_getp_material(el);
    double const zatom    = MaterialData_get_Z(material);
    double const emr      = MaterialData_get_nuclear_radius(material);
    double const hcut     = MaterialData_get_hcut(material);

    EverestRandomData evran = EverestCollimatorData_getp_random_generator(el);
    EverestRandomData_set_rutherford(evran, zatom, emr, hcut);

    //start_per_particle_block (part0->part)
        int8_t is_tracking = xcoll_assert_tracking(part);
        int8_t rng_set     = xcoll_assert_rng_set(part);

        if ( is_tracking && rng_set ) {
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
        drift_6d(part0, inactive_front+active_length+inactive_back);
    } else {
        drift_6d(part0, inactive_front);
        track_collimator(el, part0);
        //TODO: check that this only drifts surviving particles!!!
        drift_6d(part0, inactive_back);
    }
}

#endif
