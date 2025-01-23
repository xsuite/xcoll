// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_H
#define XCOLL_EVEREST_CRYSTAL_H
#include <math.h>
#include <stdio.h>


/*gpufun*/
int8_t EverestCrystalData_get_record_impacts(EverestCrystalData el){
    return EverestCrystalData_get__record_interactions(el) % 2;
}

/*gpufun*/
int8_t EverestCrystalData_get_record_exits(EverestCrystalData el){
    return (EverestCrystalData_get__record_interactions(el) >> 1) % 2;
}

/*gpufun*/
int8_t EverestCrystalData_get_record_scatterings(EverestCrystalData el){
    return (EverestCrystalData_get__record_interactions(el) >> 2) % 2;
}


void EverestCrystal_set_material(EverestCrystalData el){
    CrystalMaterialData material = EverestCrystalData_getp__material(el);
    RandomRutherfordData rng = EverestCrystalData_getp_rutherford_rng(el);
    RandomRutherford_set_by_xcoll_material(rng, (GeneralMaterialData) material);
}


/*gpufun*/
CrystalGeometry EverestCrystal_init_geometry(EverestCrystalData el, LocalParticle* part0){
    CrystalGeometry cg = (CrystalGeometry) malloc(sizeof(CrystalGeometry_));
    cg->length = EverestCrystalData_get_length(el);
    cg->side   = EverestCrystalData_get__side(el);
    if (cg->side == 0){
        kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
        return cg;
    }
    double R   = EverestCrystalData_get__bending_radius(el);
    double t_R = EverestCrystalData_get__bending_angle(el);
    cg->bending_radius = R;
    cg->bending_angle  = t_R;
    cg->miscut_angle   = EverestCrystalData_get_miscut(el);
    cg->width  = EverestCrystalData_get__width(el);
    cg->height = EverestCrystalData_get__height(el);
    cg->jaw_U  = EverestCrystalData_get__jaw_U(el);
    cg->sin_z  = EverestCrystalData_get__sin_z(el);
    cg->cos_z  = EverestCrystalData_get__cos_z(el);
    cg->sin_y  = EverestCrystalData_get__sin_y(el);
    cg->cos_y  = EverestCrystalData_get__cos_y(el);
    // Segments
    if (cg->side == 1){
        cg->segments = create_crystal(cg->bending_radius, cg->width, cg->length, cg->jaw_U, \
                                        cg->sin_y, cg->cos_y);
    } else if (cg->side == -1){
        // jaw_U is the inner corner (shifted if right-sided crystal)
        cg->segments = create_crystal(cg->bending_radius, cg->width, cg->length, cg->jaw_U - cg->width, \
                                        cg->sin_y, cg->cos_y);
    }
    // // Jaw frame is always left-sided
    // cg->segments_jf = create_crystal(cg->bending_radius, cg->width, cg->length, 0, 0, 1);
    // Bend centre
    cg->s_B = 0;
    cg->x_B = cg->bending_radius;
    // Miscut centre
    cg->s_P = -R*sin(cg->miscut_angle);
    cg->x_P = R*cos(cg->miscut_angle);
    if (cg->side == 1 && R < 0){
        // If R<0, a left-sided crystal bends towards the beam
        cg->x_P = cg->x_P + cg->width;
        cg->x_B = cg->x_B + cg->width;
    } else if (cg->side == -1 && R > 0){
        // If R>0, a right-sided crystal bends towards the beam
        cg->x_P = cg->x_P - cg->width;
        cg->x_B = cg->x_B - cg->width;
    }
    if (cg->side == -1){
        // Mirror the crystal geometry
        cg->bending_radius = -cg->bending_radius;
        cg->bending_angle  = -cg->bending_angle;
        cg->miscut_angle   = -cg->miscut_angle;
        cg->x_P            = -cg->x_P;
        cg->x_B            = -cg->x_B;
    }
    // From here on, crystal geometry parameters can always be treated as left-sided.
    // Note that the segments are not mirrored, which is fine as get_s_of_first_crossing_with_vlimit
    // is absolute (not in the jaw reference frame). It is only after a hit is registered, that we
    // need to transform the particle to the jaw reference frame.
    double Rb;
    if (cg->miscut_angle > 0){
        Rb = R - cg->width;
    } else {
        Rb = R;
    }
    cg->t_VImax = atan( (Rb*sin(t_R) - cg->s_P) / (R - Rb*cos(t_R) - cg->x_P) );
    // Impact table
    cg->record = EverestCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record){
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = EverestCrystalData_get_record_impacts(el);
        cg->record_exits = EverestCrystalData_get_record_exits(el);
    }
    return cg;
}

/*gpufun*/
void EverestCrystal_free(CrystalGeometry restrict cg){
    destroy_crystal(cg->segments);
    free(cg);
}


// TODO: it would be great if we could set EverestData as an xofield, because then we could
// run this function at creation of the collimator instead of every turn
/*gpufun*/
EverestCollData EverestCrystal_init(EverestCrystalData el, LocalParticle* part0){
    EverestCollData coll = (EverestCollData) malloc(sizeof(EverestCollData_));
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
    coll->only_mcs = CrystalMaterialData_get__only_mcs(material);
    coll->dlri     = CrystalMaterialData_get_crystal_radiation_length(material);
    coll->dlyi     = CrystalMaterialData_get_crystal_nuclear_length(material);
    coll->ai       = CrystalMaterialData_get_crystal_plane_distance(material);
    coll->eum      = CrystalMaterialData_get_crystal_potential(material);
    coll->collnt   = CrystalMaterialData_get_nuclear_collision_length(material);
    coll->eta      = 0.9;  // Hard-coded channeling saturation factor
    coll->orient   = EverestCrystalData_get__orient(el);
    // Impact table
    coll->record = EverestCrystalData_getp_internal_record(el, part0);
    coll->record_index = NULL;
    if (coll->record){
        coll->record_index = InteractionRecordData_getp__index(coll->record);
        coll->record_scatterings = EverestCrystalData_get_record_scatterings(el);
    }
    return coll;
}


/*gpufun*/
EverestData EverestCrystal_init_data(LocalParticle* part, EverestCollData restrict coll, CrystalGeometry restrict cg){
    EverestData everest = (EverestData) malloc(sizeof(EverestData_));
    everest->coll = coll;
    everest->rescale_scattering = 1;
    // Preinitialise scattering parameters
    double energy = LocalParticle_get_energy(part) / 1e9; // energy in GeV
    calculate_scattering(everest, energy);
    calculate_ionisation_properties(everest, energy);
    calculate_critical_angle(everest, part, cg, energy);
    calculate_VI_parameters(everest, part, energy);
    return everest;
}


/*gpufun*/
void EverestCrystal_track_local_particle(EverestCrystalData el, LocalParticle* part0) {
    int8_t active = EverestCrystalData_get_active(el);
    active       *= EverestCrystalData_get__tracking(el);
    double length = EverestCrystalData_get_length(el);

    // Initialise collimator data
    EverestCollData coll;
    CrystalGeometry cg;
    if (active){
        // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
        coll = EverestCrystal_init(el, part0);
        cg   = EverestCrystal_init_geometry(el, part0);

        // For info
        double const e0 = LocalParticle_get_energy0(part0);
        double t_c0  = _critical_angle0(coll, e0);
        double Rcrit = _critical_radius(coll, e0);
        double t_c = _critical_angle(coll, t_c0, Rcrit / fabs(cg->bending_radius));
        EverestCrystalData_set__critical_radius(el, Rcrit);
        EverestCrystalData_set__critical_angle(el, t_c);

        if (cg->width==0 || cg->height==0 || cg->bending_radius==0){
            kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
        }
    }

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = xcoll_check_particle_init(coll->rng, part);

            if (is_valid) {
                // Store s-location of start of crystal
                double const s_coll = LocalParticle_get_s(part);
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

                // Check if hit on jaws
                int8_t is_hit = hit_crystal_check_and_transform(part, cg);

                if (is_hit != 0) {
                    // Hit one of the jaws, so scatter
                    double remaining_length = length - LocalParticle_get_s(part);
                    // Scatter
                    EverestData everest = EverestCrystal_init_data(part, coll, cg);
                    calculate_initial_angle(everest, part, cg);
#ifdef XCOLL_USE_EXACT
                    double const xp = LocalParticle_get_exact_xp(part);
#else
                    double const xp = LocalParticle_get_xp(part);
#endif
                    if (fabs(xp - everest->t_I) < everest->t_c) {
                        energy_out = Channel(everest, part, cg, energy_in/1.e9, remaining_length)*1.e9;
                    } else {
                        energy_out = Amorphous(everest, part, cg, energy_in/1.e9, remaining_length, 1)*1.e9;
                    }
                    free(everest);
                }

                // Transform back to the lab frame
                hit_crystal_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);

                LocalParticle_set_zeta(part, zeta_in);
                // Hit and survived particles need correcting:
                if (is_hit!=0 && LocalParticle_get_state(part)>0){
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
    if (active){
        EverestCrystal_free(cg);
        free(coll);
    }
}


#endif /* XCOLL_EVEREST_CRYSTAL_H */
