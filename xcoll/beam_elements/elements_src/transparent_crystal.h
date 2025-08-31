// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_TRANSPARENT_CRY_H
#define XCOLL_TRANSPARENT_CRY_H

#include <headers/track.h>
#include <headers/checks.h>
#include <xcoll/headers/particle_states.h>


GPUFUN
int8_t TransparentCrystalData_get_record_impacts(TransparentCrystalData el){
    return TransparentCrystalData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t TransparentCrystalData_get_record_exits(TransparentCrystalData el){
    return (TransparentCrystalData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t TransparentCrystalData_get_record_scatterings(TransparentCrystalData el){
    return (TransparentCrystalData_get__record_interactions(el) >> 2) % 2;
}


GPUFUN
CrystalGeometry TransparentCrystal_init_geometry(TransparentCrystalData el, LocalParticle* part0){
    CrystalGeometry cg = (CrystalGeometry) malloc(sizeof(CrystalGeometry_));
    cg->length = TransparentCrystalData_get_length(el);
    cg->side   = TransparentCrystalData_get__side(el);
    cg->bending_radius = TransparentCrystalData_get__bending_radius(el);
    cg->bending_angle  = TransparentCrystalData_get__bending_angle(el);
    cg->width  = TransparentCrystalData_get__width(el);
    cg->height = TransparentCrystalData_get__height(el);
    cg->jaw_U  = TransparentCrystalData_get__jaw_U(el);
    cg->sin_z  = TransparentCrystalData_get__sin_z(el);
    cg->cos_z  = TransparentCrystalData_get__cos_z(el);
    cg->sin_y  = TransparentCrystalData_get__sin_y(el);
    cg->cos_y  = TransparentCrystalData_get__cos_y(el);
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
    cg->record = TransparentCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record){
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = TransparentCrystalData_get_record_impacts(el);
        cg->record_exits = TransparentCrystalData_get_record_exits(el);
    }
    // Not needed, set to zero
    cg->miscut_angle = 0;
    cg->s_P = 0;
    cg->x_P = 0;
    cg->t_VImax = 0;
    return cg;
}

GPUFUN
void TransparentCrystal_free(CrystalGeometry restrict cg){
    destroy_crystal(cg->segments);
    free(cg);
}


GPUFUN
void TransparentCrystal_track_local_particle(TransparentCrystalData el, LocalParticle* part0){
    int8_t active = TransparentCrystalData_get_active(el);
    active       *= TransparentCrystalData_get__tracking(el);
    double const length = TransparentCrystalData_get_length(el);

    // Get geometry
    CrystalGeometry cg;
    if (active){
        cg = TransparentCrystal_init_geometry(el, part0);
        if (cg->width==0 || cg->height==0 || cg->bending_radius==0){
            kill_all_particles(part0, XC_ERR_INVALID_XOFIELD);
        }
    }

    START_PER_PARTICLE_BLOCK(part0, part);
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

                // Transform back to the lab frame
                hit_crystal_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);
            }
        }
    END_PER_PARTICLE_BLOCK;
    if (active){
        TransparentCrystal_free(cg);
    }
}

#endif /* XCOLL_TRANSPARENT_CRY_H */
