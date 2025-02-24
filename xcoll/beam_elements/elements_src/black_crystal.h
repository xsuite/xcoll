// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_ABSORBER_CRY_H
#define XCOLL_ABSORBER_CRY_H


/*gpufun*/
int8_t BlackCrystalData_get_record_impacts(BlackCrystalData el){
    return BlackCrystalData_get__record_interactions(el) % 2;
}

/*gpufun*/
int8_t BlackCrystalData_get_record_exits(BlackCrystalData el){
    return (BlackCrystalData_get__record_interactions(el) >> 1) % 2;
}

/*gpufun*/
int8_t BlackCrystalData_get_record_scatterings(BlackCrystalData el){
    return (BlackCrystalData_get__record_interactions(el) >> 2) % 2;
}


/*gpufun*/
CrystalGeometry BlackCrystal_init_geometry(BlackCrystalData el, LocalParticle* part0){
    CrystalGeometry cg = (CrystalGeometry) malloc(sizeof(CrystalGeometry_));
    cg->length = BlackCrystalData_get_length(el);
    cg->side   = BlackCrystalData_get__side(el);
    cg->bending_radius = BlackCrystalData_get__bending_radius(el);
    cg->bending_angle  = BlackCrystalData_get__bending_angle(el);
    cg->width  = BlackCrystalData_get__width(el);
    cg->height = BlackCrystalData_get__height(el);
    cg->jaw_U  = BlackCrystalData_get__jaw_U(el);
    cg->sin_z  = BlackCrystalData_get__sin_z(el);
    cg->cos_z  = BlackCrystalData_get__cos_z(el);
    cg->sin_y  = BlackCrystalData_get__sin_y(el);
    cg->cos_y  = BlackCrystalData_get__cos_y(el);
    double jaw;
    if (cg->side == 1){
        jaw = cg->jaw_U;
    } else if (cg->side == -1){
        jaw = cg->jaw_U - cg->width;   // To ensure that jaw_U is the inner corner
    } else {
        kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
        return cg;
    }
    cg->segments = create_crystal(cg->bending_radius, cg->width, cg->length, jaw, cg->sin_y, cg->cos_y);
    // Impact table
    cg->record = BlackCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record){
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = BlackCrystalData_get_record_impacts(el);
        cg->record_exits = BlackCrystalData_get_record_exits(el);
    }
    // Not needed, set to zero
    cg->miscut_angle = 0;
    cg->s_P = 0;
    cg->x_P = 0;
    cg->t_VImax = 0;
    return cg;
}

/*gpufun*/
void BlackCrystal_free(CrystalGeometry restrict cg){
    destroy_crystal(cg->segments);
    free(cg);
}


/*gpufun*/
void BlackCrystal_track_local_particle(BlackCrystalData el, LocalParticle* part0){
    int8_t active = BlackCrystalData_get_active(el);
    active       *= BlackCrystalData_get__tracking(el);
    double const length = BlackCrystalData_get_length(el);

    // Get geometry
    CrystalGeometry cg;
    int8_t record_scatterings;
    if (active){
        cg = BlackCrystal_init_geometry(el, part0);
        record_scatterings = BlackCrystalData_get_record_scatterings(el);

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
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_tracking) {
                // Store s-location of start of collimator
                double s_coll = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Check if hit on jaws
                int8_t is_hit = hit_crystal_check_and_transform(part, cg);

                if (is_hit != 0){
                    LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);
                    if (record_scatterings) {
                        InteractionRecordData_log(cg->record, cg->record_index, part, XC_ABSORBED);  // In coll jaw reference frame
                    }
                }

                // Transform back to the lab frame
                hit_crystal_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);
            }
        }
    //end_per_particle_block
    if (active){
        BlackCrystal_free(cg);
    }
}

#endif /* XCOLL_ABSORBER_CRY_H */
