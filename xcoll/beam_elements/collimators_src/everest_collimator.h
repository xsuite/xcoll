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

    double tilt_L = asin(EverestCollimatorData_get__sin_yL(el));
    double tilt_R = asin(EverestCollimatorData_get__sin_yR(el));
    double sin_zL = asin(EverestCollimatorData_get__sin_zL(el));
    double cos_zL = asin(EverestCollimatorData_get__cos_zL(el));
    double jaw_LU = EverestCollimatorData_get__jaw_LU(el);
    double jaw_RU = EverestCollimatorData_get__jaw_RU(el);
    double jaw_LD = EverestCollimatorData_get__jaw_LD(el);
    double jaw_RD = EverestCollimatorData_get__jaw_RD(el);
    double jaw_L = (jaw_LU + jaw_LD)/2;
    double jaw_R = (jaw_RU + jaw_RD)/2;

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                // Save the end position of the particle
                double s_coll = LocalParticle_get_s(part);

                // Store initial coordinates for updating later
                double const rpp_in  = LocalParticle_get_rpp(part);
                double const rvv_in  = LocalParticle_get_rvv(part);
                double const e0      = LocalParticle_get_energy0(part) / 1e9; // Reference energy in GeV
                double const beta0   = LocalParticle_get_beta0(part);
                double const ptau_in = LocalParticle_get_ptau(part);
                double const x_in    = LocalParticle_get_x(part);
                double const px      = LocalParticle_get_px(part);
                double const y_in    = LocalParticle_get_y(part);
                double const py      = LocalParticle_get_py(part);
                double const zeta_in = LocalParticle_get_zeta(part);
                double const px_in   =  cos_zL * px + sin_zL * py;
                double const py_in   = -sin_zL * px + cos_zL * py;
                double p0 = LocalParticle_get_p0c(part) / 1e9;
                // TODO: missing correction due to m/m0 (but also wrong in xpart...)
                double energy = p0*ptau_in + e0; // energy, not momentum, in GeV
                int is_abs = 0;

                // Check if hit on jaws
                int8_t is_hit = hit_jaws_check_and_transform(part, cg);

                // Hit one of the jaws, so scatter
                double s_part = LocalParticle_get_s(part) - s_coll;
                double remaining_length = length - s_part;

                if (is_hit == 1){
                    // left jaw
                    XYShift_single_particle(part, jaw_L, 0);
                     // Include collimator tilt
                    double rot_shift = YRotation_single_particle_rotate_only(part, s_part, tilt_L);
                    if (fabs(s_part) < 1.e-10){
                        Drift_single_particle_4d(part, -rot_shift);
                    }
                } else if (is_hit == -1){
                    // Right jaw
                    XYShift_single_particle(part, jaw_R, 0);
                    LocalParticle_scale_x(part, -1);
                    LocalParticle_scale_px(part, -1);
                     // Include collimator tilt
                    double rot_shift = YRotation_single_particle_rotate_only(part, s_part, -tilt_R);
                    if (fabs(s_part) < 1.e-10){
                        Drift_single_particle_4d(part, -rot_shift);
                    }
                }

                if (is_hit != 0) {
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
                        Drift_single_particle_4d(part, remaining_length-s_out);
                    }
                }


                if (is_hit == 1){
                    double rot_shift = YRotation_single_particle_rotate_only(part, length, -tilt_L); // Wrong, should be final s-position (in case particle is absorbed..)
                    Drift_single_particle_4d(part, length-rot_shift);
                    XYShift_single_particle(part, -jaw_LU, 0);
                } else if (is_hit == -1){
                    double rot_shift = YRotation_single_particle_rotate_only(part, length, tilt_R);
                    Drift_single_particle_4d(part, length-rot_shift);
                    LocalParticle_scale_x(part, -1);
                    LocalParticle_scale_px(part, -1);
                    XYShift_single_particle(part, -jaw_RU, 0);
                }

                // Update energy    ---------------------------------------------------
                // Only particles that hit the jaw and survived need to be updated
                if (is_hit!=0 && is_abs==0){
                    double ptau_out = (energy - e0) / (e0 * beta0);
                    LocalParticle_update_ptau(part, ptau_out);
                }
            
                // Update 4D coordinates    -------------------------------------------
                // Absorbed particles get their coordinates set to the entrance of collimator
                if (is_abs>0){
                    LocalParticle_set_x(part, x_in);
                    LocalParticle_set_px(part, px_in);
                    LocalParticle_set_y(part, y_in);
                    LocalParticle_set_py(part, py_in);
                }

                LocalParticle_set_zeta(part, zeta_in);
                // Hit and survived particles need correcting:
                if (is_hit!=0 && is_abs==0){
                    double px  = LocalParticle_get_px(part);
                    double py  = LocalParticle_get_py(part);
                    double rvv = LocalParticle_get_rvv(part);
                    double rpp = LocalParticle_get_rpp(part);
                    // First we drift half the length with the old angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, px_in*rpp_in, py_in*rpp_in, length/2) );
                    // then half the length with the new angles:
                    LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, px*rpp, py*rpp, length/2) );
                }
            
                // Update state    ----------------------------------------------------
                if (is_abs > 0){
                    LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
                }

                // Transform back to the lab frame
                hit_jaws_transform_back(is_hit, part, cg);
            }
        }
    //end_per_particle_block
    free(cg);
    free(coll);
}


#endif /* XCOLL_EVEREST_COLL_H */
