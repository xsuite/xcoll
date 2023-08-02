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


// TODO: it would be great if we could set EverestData as an xofield, because then we could
// run this function at creation of the collimator instead of every turn
/*gpufun*/
EverestData EverestCrystal_init(EverestCrystalData el, LocalParticle* part0){

    EverestData coll = (EverestData) malloc(sizeof(EverestData_));

    // Random generator and material
    coll->rng = EverestCrystalData_getp_rutherford_rng(el);
    CrystalMaterialData material = EverestCrystalData_getp__material(el);
    coll->exenergy = CrystalMaterialData_get_excitation_energy(material)*1.0e3; // MeV
    coll->rho      = CrystalMaterialData_get_density(material);
    coll->anuc     = CrystalMaterialData_get_A(material);
    coll->zatom    = CrystalMaterialData_get_Z(material);
    coll->bnref    = CrystalMaterialData_get_nuclear_elastic_slope(material);
    coll->csref[0] = CrystalMaterialData_get_cross_section(material, 0);
    coll->csref[1] = CrystalMaterialData_get_cross_section(material, 1);
    coll->csref[5] = CrystalMaterialData_get_cross_section(material, 5);
    coll->dlri     = CrystalMaterialData_get_crystal_radiation_length(material);
    coll->dlyi     = CrystalMaterialData_get_crystal_nuclear_length(material);
    coll->ai       = CrystalMaterialData_get_crystal_plane_distance(material);
    coll->eum      = CrystalMaterialData_get_crystal_potential(material);
    coll->collnt   = CrystalMaterialData_get_nuclear_collision_length(material);

    // Impact table
    coll->record = EverestCrystalData_getp_internal_record(el, part0);
    coll->record_index = NULL;
    if (coll->record){
        coll->record_index = CollimatorImpactsData_getp__index(coll->record);
    }

    // Geometry
    // TODO: this should in principle not be in this struct
    coll->aperture = EverestCrystalData_get_jaw_L(el) - EverestCrystalData_get_jaw_R(el);
    coll->offset   = ( EverestCrystalData_get_jaw_L(el) + EverestCrystalData_get_jaw_R(el) ) /2;
    coll->tilt_L   = asin(EverestCrystalData_get_sin_yL(el));
    coll->tilt_R   = asin(EverestCrystalData_get_sin_yR(el));
    if (fabs(coll->tilt_R) > 1.e-10){
        kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
    };
    coll->side     = EverestCrystalData_get__side(el);  // TODO: so far only left-sided crystals
    if (coll->side != 1){
        kill_all_particles(part0, XC_ERR_NOT_IMPLEMENTED);
    };
    // TODO: this should stay here
    coll->bend_r     = EverestCrystalData_get__bending_radius(el);
    coll->bend_ang   = EverestCrystalData_get__bending_angle(el);
    coll->tilt       = EverestCrystalData_get_align_angle(el) + coll->tilt_L;   // TODO: only left-sided crystals
    // double const cry_bend   = length/cry_rcurv; //final value (with corrected length)
    // THIS IS WRONG! Was a mistranslation from SixTrack 4 to SixTrack 5
    // Difference is so small that this was never caught.
    // Furthermore, we removed the adaptation of the scatter length, because
    // 1) it was implemented wrong (passed unnoticed due to small effect)
    // 2) we should not use the adapted scatter length, as we rotate the S-X frame, so
    //    we anyway have to drift the full length!    
    coll->amorphous_layer = EverestCrystalData_get_thick(el);
    coll->xdim            = EverestCrystalData_get_xdim(el);
    coll->ydim            = EverestCrystalData_get_ydim(el);
    coll->orient          = EverestCrystalData_get__orient(el);
    coll->miscut          = EverestCrystalData_get_miscut(el);

    return coll;
}


/*gpufun*/
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t active      = EverestCrystalData_get_active(el);
    active            *= EverestCrystalData_get__tracking(el);
    double const inactive_front = EverestCrystalData_get_inactive_front(el);
    double const active_length  = EverestCrystalData_get_active_length(el);
    double const inactive_back  = EverestCrystalData_get_inactive_back(el);

    // Collimator geometry
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

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    EverestData coll = EverestCrystal_init(el, part0);

    // Preinitialise scattering parameters
    double const energy0 = LocalParticle_get_energy0(&part0[0]) / 1e9; // Reference energy in GeV
    calculate_scattering(coll, energy0, 1.);
    calculate_ionisation_properties(coll, energy0);

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, inactive_front + active_length + inactive_back);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
            int8_t rng_is_set  = assert_rng_set(part, RNG_ERR_SEEDS_NOT_SET);
            int8_t ruth_is_set = assert_rutherford_set(coll->rng, part, RNG_ERR_RUTH_NOT_SET);

            if (is_tracking && rng_is_set && ruth_is_set) {
                // Drift inactive front
                Drift_single_particle(part, inactive_front);

                // Move to collimator frame
                XYShift_single_particle(part, co_x, co_y);
                SRotation_single_particle(part, sin_zL, cos_zL);

                scatter_cry(coll, part, active_length);

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

    free(coll);
}


#endif /* XCOLL_EVEREST_CRYSTAL_H */
