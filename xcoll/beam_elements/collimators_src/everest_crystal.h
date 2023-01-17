#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H

#include <math.h>

// TODO:
//    Do not split 4d and zeta in drifts
//    Use drift function from xtrack Drift element (call its C function)

/*gpufun*/
void cry_drift_6d(LocalParticle* part0, double length) {
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
void track_crystal(EverestCrystalData el, LocalParticle* part0) {

    double const energy0 = LocalParticle_get_energy0(&part0[0]) / 1e9; // Reference energy in GeV

    // Material properties
    CrystalMaterialData material = EverestCrystalData_getp_material(el);
    double const zatom    = CrystalMaterialData_get_Z(material);
    double const emr      = CrystalMaterialData_get_nuclear_radius(material);
    double const hcut     = CrystalMaterialData_get_hcut(material);

    // Calculate scattering parameters
    set_rutherford_parameters(zatom, emr, hcut);

    //start_per_particle_block (part0->part)
        scatter_cry(el, part);
    //end_per_particle_block
}

/*gpufun*/
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t const is_active      = EverestCrystalData_get__active(el);
    double const inactive_front = EverestCrystalData_get_inactive_front(el);
    double const active_length  = EverestCrystalData_get_active_length(el);
    double const inactive_back  = EverestCrystalData_get_inactive_back(el);

    if (!is_active){
        // Drift full length
        cry_drift_6d(part0, inactive_front+active_length+inactive_back);
    } else {
        cry_drift_6d(part0, inactive_front);
        track_crystal(el, part0);
        cry_drift_6d(part0, inactive_back);
    }
}

#endif
