// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_BLOCK_H
#define XCOLL_EVEREST_BLOCK_H
#include <math.h>
#include <stdio.h>


/*gpufun*/
int8_t EverestBlockData_get_record_impacts(EverestBlockData el){
    return EverestBlockData_get__record_interactions(el) % 2;
}

/*gpufun*/
int8_t EverestBlockData_get_record_exits(EverestBlockData el){
    return (EverestBlockData_get__record_interactions(el) >> 1) % 2;
}

/*gpufun*/
int8_t EverestBlockData_get_record_scatterings(EverestBlockData el){
    return (EverestBlockData_get__record_interactions(el) >> 2) % 2;
}

/*gpufun*/
void EverestBlock_set_material(EverestBlockData el){
    MaterialData material = EverestBlockData_getp__material(el);
    RandomRutherfordData rng = EverestBlockData_getp_rutherford_rng(el);
    RandomRutherford_set_by_xcoll_material(rng, (GeneralMaterialData) material);
}


/*gpufun*/
EverestCollData EverestBlock_init(EverestBlockData el, LocalParticle* part0, int8_t active){
    EverestCollData coll = (EverestCollData) malloc(sizeof(EverestCollData_));
    if (active){ // This is needed in order to avoid that the initialisation is called during a twiss!
        // Random generator and material
        coll->rng = EverestBlockData_getp_rutherford_rng(el);
        MaterialData material = EverestBlockData_getp__material(el);
        coll->exenergy = MaterialData_get_excitation_energy(material)*1.0e3; // MeV
        coll->rho      = MaterialData_get_density(material);
        coll->anuc     = MaterialData_get_A(material);
        coll->zatom    = MaterialData_get_Z(material);
        coll->bnref    = MaterialData_get_nuclear_elastic_slope(material);
        coll->radl     = MaterialData_get_radiation_length(material);
        coll->csref[0] = MaterialData_get_cross_section(material, 0);
        coll->csref[1] = MaterialData_get_cross_section(material, 1);
        coll->csref[5] = MaterialData_get_cross_section(material, 5);
        coll->only_mcs = MaterialData_get__only_mcs(material);

        // Impact table
        coll->record = EverestBlockData_getp_internal_record(el, part0);
        coll->record_index = NULL;
        coll->record_scatterings = 0;
        if (coll->record){
            coll->record_index = InteractionRecordData_getp__index(coll->record);
            coll->record_scatterings = EverestBlockData_get_record_scatterings(el);
        }
    }
    return coll;
}


/*gpufun*/
EverestData EverestBlock_init_data(LocalParticle* part, EverestCollData coll){
    EverestData everest = (EverestData) malloc(sizeof(EverestData_));
    everest->coll = coll;
    everest->rescale_scattering = 1;
    // Preinitialise scattering parameters
    double energy = LocalParticle_get_energy(part) / 1e9; // energy in GeV
    calculate_scattering(everest, energy);
    calculate_ionisation_properties(everest, energy);
    return everest;
}


/*gpufun*/
void EverestBlock_track_local_particle(EverestBlockData el, LocalParticle* part0) {
    int8_t active = EverestBlockData_get_active(el);
    double const length = EverestBlockData_get_length(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    EverestCollData coll = EverestBlock_init(el, part0, active);

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                // Store s-location of start of block
                double s_block = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Store initial coordinates for updating later
                double const rvv_in     = LocalParticle_get_rvv(part);
#ifdef XCOLL_USE_EXACT
                double const xp_in      = LocalParticle_get_exact_xp(part);
                double const yp_in      = LocalParticle_get_exact_yp(part);
#else
                double const xp_in      = LocalParticle_get_xp(part);
                double const yp_in      = LocalParticle_get_yp(part);
#endif
                double const zeta_in    = LocalParticle_get_zeta(part);
                double const energy_in  = LocalParticle_get_energy(part);
                double energy_out;

                EverestData everest = EverestBlock_init_data(part, coll);
                energy_out = jaw(everest, part, energy_in, length, 0);
                free(everest);
                LocalParticle_add_to_s(part, s_block);

                LocalParticle_set_zeta(part, zeta_in);
                // Survived particles need correcting:
                if (LocalParticle_get_state(part)>0){
                    // Update energy; the last flag keeps angles constant (even valid for exact angles!)
                    LocalParticle_add_to_energy(part, energy_out - energy_in, 0);
                    // Update zeta
#ifdef XCOLL_USE_EXACT
                    double xp  = LocalParticle_get_exact_xp(part);
                    double yp  = LocalParticle_get_exact_yp(part);
#else
                    double xp  = LocalParticle_get_xp(part);
                    double yp  = LocalParticle_get_yp(part);
#endif
                    double rvv = LocalParticle_get_rvv(part);
                    // First we drift half the length with the old angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, xp_in, yp_in, length/2) );
                    // then half the length with the new angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, xp, yp, length/2) );
                }
            }
        }
    //end_per_particle_block
    free(coll);
}



#endif /* XCOLL_EVEREST_BLOCK_H */