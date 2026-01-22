// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_FLUKA_CRY_H
#define XCOLL_FLUKA_CRY_H


/*gpufun*/
int8_t FlukaCrystalData_get_record_impacts(FlukaCrystalData el){
    return FlukaCrystalData_get__record_interactions(el) % 2;
}

/*gpufun*/
int8_t FlukaCrystalData_get_record_exits(FlukaCrystalData el){
    return (FlukaCrystalData_get__record_interactions(el) >> 1) % 2;
}

/*gpufun*/
int8_t FlukaCrystalData_get_record_scatterings(FlukaCrystalData el){
    return (FlukaCrystalData_get__record_interactions(el) >> 2) % 2;
}


/*gpufun*/
CrystalGeometry FlukaCrystal_init_geometry(FlukaCrystalData el, LocalParticle* part0){
    CrystalGeometry cg = (CrystalGeometry) malloc(sizeof(CrystalGeometry_));
    cg->length = FlukaCrystalData_get_length(el);
    cg->side   = FlukaCrystalData_get__side(el);
    cg->bending_radius = FlukaCrystalData_get__bending_radius(el);
    cg->bending_angle  = FlukaCrystalData_get__bending_angle(el);
    cg->width  = FlukaCrystalData_get__width(el);
    cg->height = FlukaCrystalData_get__height(el);
    cg->jaw_U  = FlukaCrystalData_get__jaw_U(el);
    cg->sin_z  = FlukaCrystalData_get__sin_z(el);
    cg->cos_z  = FlukaCrystalData_get__cos_z(el);
    cg->sin_y  = FlukaCrystalData_get__sin_y(el);
    cg->cos_y  = FlukaCrystalData_get__cos_y(el);
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
    cg->record = FlukaCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record){
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = FlukaCrystalData_get_record_impacts(el);
        cg->record_exits = FlukaCrystalData_get_record_exits(el);
    }
    // Not needed, set to zero
    cg->miscut_angle = 0;
    cg->s_P = 0;
    cg->x_P = 0;
    cg->t_VImax = 0;
    return cg;
}

/*gpufun*/
void FlukaCrystal_free(CrystalGeometry restrict cg){
    destroy_crystal(cg->segments);
    free(cg);
}


/*gpufun*/
void FlukaCrystal_track_local_particle(FlukaCrystalData el, LocalParticle* part0){
    int8_t active = FlukaCrystalData_get_active(el);
    active       *= FlukaCrystalData_get__tracking(el);

    // Get geometry
    CrystalGeometry cg;
    if (active){
        cg = FlukaCrystal_init_geometry(el, part0);
        if (cg->width==0 || cg->height==0 || cg->bending_radius==0){
            kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
        }
    }

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            double length = FlukaCrystalData_get_length(el);
            length += FlukaCrystalData_get_length_front(el);
            length += FlukaCrystalData_get_length_back(el);
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_tracking) {
                // Move to start of crystal
                double length_front = FlukaCrystalData_get_length_front(el);
#ifdef XCOLL_USE_EXACT
                Drift_single_particle_exact(part, length_front);
#else
                Drift_single_particle_expanded(part, length_front);
#endif

                // Store s-location of start of crystal
                double s_coll = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Check if hit on jaws
                int8_t is_hit;
                if (cg->record_impacts){
                    // Check with transformation (to log impacts).
                    // Particle will be at impact position or at end.
                    is_hit = hit_crystal_check_and_transform(part, cg);
                    if (is_hit != 0){
                        // Particle needs to return to start position.
                        hit_crystal_return(is_hit, part, cg);
                    }

                } else {
                    // Only check. Particle will be at start position.
                    is_hit = hit_crystal_check(part, cg);
                    if (is_hit == 0){
                        // Drift to end.
#ifdef XCOLL_USE_EXACT
                        Drift_single_particle_exact(part, cg->length);
#else
                        Drift_single_particle_expanded(part, cg->length);
#endif
                    }
                }

                if (is_hit == 0){
                    // Drift to end.
                    double length_back = FlukaCrystalData_get_length_back(el);
#ifdef XCOLL_USE_EXACT
                    Drift_single_particle_exact(part, length_back);
#else
                    Drift_single_particle_expanded(part, length_back);
#endif

                } else {
                    // Mark for FLUKA processing.
                    LocalParticle_set_state(part, XC_HIT_ON_FLUKA_COLL);
                    // Return to start position
#ifdef XCOLL_USE_EXACT
                    Drift_single_particle_exact(part, -length_front);
#else
                    Drift_single_particle_expanded(part, -length_front);
#endif
                }
                LocalParticle_add_to_s(part, s_coll);
            }
        }
    //end_per_particle_block
    if (active){
        FlukaCrystal_free(cg);
    }
}

#endif /* XCOLL_FLUKA_CRY_H */
