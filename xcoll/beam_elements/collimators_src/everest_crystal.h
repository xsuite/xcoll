// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H
#include <math.h>
#include <stdio.h>


void EverestCrystal_set_material(EverestCrystalData el){
    CrystalMaterialData material = EverestCrystalData_getp__material(el);
    RandomRutherfordData rng = EverestCrystalData_getp_rutherford_rng(el);
    RandomRutherford_set_by_xcoll_material(rng, (GeneralMaterialData) material);
}


// TODO: it would be great if we could set EverestData as an xofield, because then we could
// run this function at creation of the collimator instead of every turn
/*gpufun*/
EverestCollData EverestCrystal_init(EverestCrystalData el, LocalParticle* part0, int8_t active){
    EverestCollData coll = (EverestCollData) malloc(sizeof(EverestCollData_));
    if (active){ // This is needed in order to avoid that the initialisation is called during a twiss!
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
        double jaw_L = (EverestCrystalData_get__jaw_LU(el) + EverestCrystalData_get__jaw_LD(el))/2.;
        double jaw_R = (EverestCrystalData_get__jaw_RU(el) + EverestCrystalData_get__jaw_RD(el))/2.;
        coll->aperture = jaw_L - jaw_R;
        coll->offset   = (jaw_L + jaw_R) /2;
        coll->tilt_L   = asin(EverestCrystalData_get__sin_yL(el));
        coll->tilt_R   = asin(EverestCrystalData_get__sin_yR(el));
        if (fabs(coll->tilt_R) > 1.e-10){
            printf("Crystals have to be left-sided for now, so tilt_R should not be set!");
            fflush(stdout);
            kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
        };
        coll->side     = EverestCrystalData_get__side(el);  // TODO: so far only left-sided crystals
        if (coll->side != 1){
            printf("Crystals have to be left-sided for now!");
            fflush(stdout);
            kill_all_particles(part0, XC_ERR_NOT_IMPLEMENTED);
        };
        // TODO: this should stay here
        double R         = EverestCrystalData_get__bending_radius(el);
        coll->bend_r     = R;
        double t_R       = EverestCrystalData_get__bending_angle(el);
        coll->bend_ang   = t_R;
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
        coll->s_P             = -coll->bend_r*sin(coll->miscut);
        coll->x_P             = coll->bend_r*cos(coll->miscut);
        double Rb;
        if (coll->miscut >0){
            Rb = R - coll->xdim;
        } else {
            Rb = R;
        }
        coll->t_VImax = atan( (Rb*sin(t_R) - coll->s_P) / (R - Rb*cos(t_R) - coll->x_P) );
    }
    return coll;
}


/*gpufun*/
EverestData EverestCrystal_init_data(LocalParticle* part, EverestCollData coll){
    EverestData everest = (EverestData) malloc(sizeof(EverestData_));
    everest->coll = coll;
    everest->rescale_scattering = 1;
#ifndef XCOLL_REFINE_ENERGY
    // Preinitialise scattering parameters
    double charge_ratio = LocalParticle_get_charge_ratio(part);
    double mass_ratio = charge_ratio / LocalParticle_get_chi(part);
    double energy = ( LocalParticle_get_ptau(part) + 1 / LocalParticle_get_beta0(part)
                     ) * mass_ratio * LocalParticle_get_p0c(part) / 1e9; // energy in GeV
    calculate_scattering(everest, energy);
    calculate_ionisation_properties(everest, energy);
    calculate_critical_angle(everest, part, energy);
    calculate_VI_parameters(everest, part, energy);
#endif
    return everest;
}


/*gpufun*/
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t active = EverestCrystalData_get_active(el);
    active       *= EverestCrystalData_get__tracking(el);
    double length = EverestCrystalData_get_length(el);

    // TODO: use xtrack C-code for rotation element
    // TODO: we are ignoring the angle of the right jaw
    // TODO: is a crystal always one-sided...?
    double const sin_zL     = EverestCrystalData_get__sin_zL(el);
    double const cos_zL     = EverestCrystalData_get__cos_zL(el);
    double const sin_zR     = EverestCrystalData_get__sin_zR(el);
    double const cos_zR     = EverestCrystalData_get__cos_zR(el);
    if (fabs(sin_zL-sin_zR) > 1.e-10 || fabs(cos_zL-cos_zR) > 1.e-10 ){
        printf("Jaws with different angles not yet implemented!");
        fflush(stdout);
        kill_all_particles(part0, XC_ERR_NOT_IMPLEMENTED);
    };

    double t_c = 0;

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    EverestCollData coll = EverestCrystal_init(el, part0, active);

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                double const s_coll = LocalParticle_get_s(part);

                // Move to collimator frame
                SRotation_single_particle(part, sin_zL, cos_zL);
                XYShift_single_particle(part, coll->offset, 0);

                EverestData everest = EverestCrystal_init_data(part0, coll);
                scatter_cry(everest, part, length);

                // Temporary workaround to store the critical angle for use later
                double energy0 = LocalParticle_get_energy0(part)/1.e9;
                calculate_critical_angle(everest, part, energy0);
                t_c = everest->t_c;
                free(everest);

                // Return from collimator frame
                XYShift_single_particle(part, -coll->offset, 0);
                SRotation_single_particle(part, -sin_zL, cos_zL);

                // Surviving particles are put at same numerical s, and drifted inactive back
                if (LocalParticle_get_state(part) > 0){
                    LocalParticle_set_s(part, s_coll + length);
                }
            }
        }
    //end_per_particle_block

    free(coll);
    EverestCrystalData_set__critical_angle(el, t_c);
}


#endif /* XCOLL_EVEREST_CRYSTAL_H */
