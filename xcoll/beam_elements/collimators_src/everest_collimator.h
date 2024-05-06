// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_COLL_H
#define XCOLL_EVEREST_COLL_H
#include <math.h>
#include <stdio.h>


void EverestCollimator_set_material(EverestCollimatorData el){
    MaterialData material = EverestCollimatorData_getp__material(el);
    RandomRutherfordData rng = EverestCollimatorData_getp_rutherford_rng(el);
    RandomRutherford_set_by_xcoll_material(rng, (GeneralMaterialData) material);
}


/*gpufun*/
CollimatorGeometry EverestCollimator_init_geometry(EverestCollimatorData el, LocalParticle* part0, int8_t active){
    CollimatorGeometry cg = (CollimatorGeometry) malloc(sizeof(CollimatorGeometry_));
    if (active){ // This is needed in order to avoid that the initialisation is called during a twiss!
        // Collimator jaws (with tilts)
        cg->jaw_LU = EverestCollimatorData_get__jaw_LU(el);
        cg->jaw_LD = EverestCollimatorData_get__jaw_LD(el);
        cg->jaw_RU = EverestCollimatorData_get__jaw_RU(el);
        cg->jaw_RD = EverestCollimatorData_get__jaw_RD(el);
        // TODO: need shortening of active length!
        cg->length = EverestCollimatorData_get_length(el);
        cg->side   = EverestCollimatorData_get__side(el);
        // Get angles of jaws
        cg->sin_zL = EverestCollimatorData_get__sin_zL(el);
        cg->cos_zL = EverestCollimatorData_get__cos_zL(el);
        cg->sin_zR = EverestCollimatorData_get__sin_zR(el);
        cg->cos_zR = EverestCollimatorData_get__cos_zR(el);
        cg->sin_zDiff = EverestCollimatorData_get__sin_zDiff(el);
        cg->cos_zDiff = EverestCollimatorData_get__cos_zDiff(el);
        cg->jaws_parallel = EverestCollimatorData_get__jaws_parallel(el);
        // Tilts
        cg->sin_yL = EverestCollimatorData_get__sin_yL(el);
        cg->cos_yL = EverestCollimatorData_get__cos_yL(el);
        cg->tan_yL = EverestCollimatorData_get__tan_yL(el);
        cg->sin_yR = EverestCollimatorData_get__sin_yR(el);
        cg->cos_yR = EverestCollimatorData_get__cos_yR(el);
        cg->tan_yR = EverestCollimatorData_get__tan_yR(el);
        // Impact table:  need it here to record touches
        cg->record = EverestCollimatorData_getp_internal_record(el, part0);
        cg->record_index = NULL;
        cg->record_touches = 0;
        if (cg->record){
            cg->record_index = InteractionRecordData_getp__index(cg->record);
            cg->record_touches = EverestCollimatorData_get_record_touches(el);
        }
    }

    return cg;
}


// TODO: it would be great if we could set EverestData as an xofield, because then we could
// run this function at creation of the collimator instead of every turn
// Hmmmm this should be called whenever we change an xofield
/*gpufun*/
EverestCollData EverestCollimator_init(EverestCollimatorData el, LocalParticle* part0, int8_t active){
    EverestCollData coll = (EverestCollData) malloc(sizeof(EverestCollData_));
    if (active){ // This is needed in order to avoid that the initialisation is called during a twiss!
        // Random generator and material
        coll->rng = EverestCollimatorData_getp_rutherford_rng(el);
        MaterialData material = EverestCollimatorData_getp__material(el);
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
        // Impact table:  need it here to record interactions
        coll->record = EverestCollimatorData_getp_internal_record(el, part0);
        coll->record_index = NULL;
        coll->record_scatterings = 0;
        if (coll->record){
            coll->record_index = InteractionRecordData_getp__index(coll->record);
            coll->record_scatterings = EverestCollimatorData_get_record_scatterings(el);
        }
    }

    return coll;
}


/*gpufun*/
EverestData EverestCollimator_init_data(LocalParticle* part, EverestCollData coll){
    EverestData everest = (EverestData) malloc(sizeof(EverestData_));
    everest->coll = coll;
    everest->rescale_scattering = 1;
#ifndef XCOLL_REFINE_ENERGY
    // Preinitialise scattering parameters
    double charge_ratio = LocalParticle_get_charge_ratio(part);
    double mass_ratio = charge_ratio / LocalParticle_get_chi(part);
    double energy = ( LocalParticle_get_ptau(part) + 1 / LocalParticle_get_beta0(part)
                     ) * mass_ratio * LocalParticle_get_p0c(part) / 1e9; // energy in GeV
    energy = LocalParticle_get_energy0(part) / 1e9;
    calculate_scattering(everest, energy);
    calculate_ionisation_properties(everest, energy);
#endif
    return everest;
}


/*gpufun*/
void EverestCollimator_track_local_particle(EverestCollimatorData el, LocalParticle* part0) {
    int8_t active = EverestCollimatorData_get_active(el);
    active       *= EverestCollimatorData_get__tracking(el);
    double const length = EverestCollimatorData_get_length(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    EverestCollData coll  = EverestCollimator_init(el, part0, active);
    CollimatorGeometry cg = EverestCollimator_init_geometry(el, part0, active);

    double sin_zL = asin(EverestCollimatorData_get__sin_zL(el));
    double cos_zL = asin(EverestCollimatorData_get__cos_zL(el));

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                // Store s-location of start of collimator
                double s_coll = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Store initial coordinates for updating later
                double const rpp_in  = LocalParticle_get_rpp(part);
                double const rvv_in  = LocalParticle_get_rvv(part);
                double const e0      = LocalParticle_get_energy0(part) / 1e9; // Reference energy in GeV
                double const beta0   = LocalParticle_get_beta0(part);
                double const ptau_in = LocalParticle_get_ptau(part);
                double const px      = LocalParticle_get_px(part);
                double const py      = LocalParticle_get_py(part);
                double const zeta_in = LocalParticle_get_zeta(part);
                double const xp_in   = (cos_zL * px + sin_zL * py)*rpp_in;
                double const yp_in   = (-sin_zL * px + cos_zL * py)*rpp_in;
                double p0 = LocalParticle_get_p0c(part) / 1e9;
                // TODO: missing correction due to m/m0 (but also wrong in xpart...)
                double energy = p0*ptau_in + e0; // energy, not momentum, in GeV
                int is_abs = 0;

                // Check if hit on jaws
                int8_t is_hit = hit_jaws_check_and_transform(part, cg);

                if (is_hit != 0) {
                    // Hit one of the jaws, so scatter
                    double remaining_length = length - LocalParticle_get_s(part);
                    // Scatter
                    EverestData everest = EverestCollimator_init_data(part, coll);
                    double* jaw_result = jaw(everest, part, energy, remaining_length, 1);
                    free(everest);
                    energy = jaw_result[0];
                    if (jaw_result[1] == 1){
                        is_abs = 1;
                    }
                    double s_out = jaw_result[2];
                    free(jaw_result);
        
                    if (is_abs != 1) {
                       // Do the rest drift, if particle left collimator early
                        Drift_single_particle(part, remaining_length-s_out);
                    }
                }
            
                // Update state    ----------------------------------------------------
                if (is_abs > 0){
                    LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
                }

                // Transform back to the lab frame
                hit_jaws_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);

                // Update energy    ---------------------------------------------------
                // Only particles that hit the jaw and survived need to be updated
                if (is_hit!=0 && is_abs==0){
                    double ptau_out = (energy - e0) / (e0 * beta0);
                    LocalParticle_update_ptau(part, ptau_out);
                }

                LocalParticle_set_zeta(part, zeta_in);
                // Hit and survived particles need correcting:
                if (is_hit!=0 && is_abs==0){
                    double px  = LocalParticle_get_px(part);
                    double py  = LocalParticle_get_py(part);
                    double rvv = LocalParticle_get_rvv(part);
                    double rpp = LocalParticle_get_rpp(part);
                    // First we drift half the length with the old angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, xp_in, yp_in, length/2) );
                    // then half the length with the new angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, px*rpp, py*rpp, length/2) );
                }
            }
        }
    //end_per_particle_block
    free(cg);
    free(coll);
}


#endif /* XCOLL_EVEREST_COLL_H */
