// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_CRY_GEOM_H
#define XCOLL_CRY_GEOM_H
#include <math.h>
#include <stdio.h>

typedef struct CrystalGeometry_ {
    // Crystal inner upstream corner (with tilt)
    double jaw_U;
    double length;
    int8_t side;
    // Jaw angle
    double sin_z;
    double cos_z;
    // Tilt
    double sin_y;
    double cos_y;
    // Crystal geometry
    double bending_radius;
    double bending_angle; // probs not needed
    double width;
    double height;
    // Segments
    Segment* segments;
    // Impact table
    InteractionRecordData record;
    RecordIndex record_index;
    int8_t record_touches;
} CrystalGeometry_;
typedef CrystalGeometry_* CrystalGeometry;


// This function checks if a particle hits a jaw (and which).
// Return value: 0 (no hit), 1 (hit on left jaw), -1 (hit on right jaw).
// Furthermore, the particle is moved to the location where it hits the jaw (drifted to the end if no hit),
//              and transformed to the reference frame of that jaw.
/*gpufun*/
int8_t hit_crystal_check_and_transform(LocalParticle* part, CrystalGeometry cg){
    double part_x, part_tan_x, part_y, part_tan_y;
    int8_t is_hit = 0;
    double s = 1.e21;

    // Crystal should be single-sided
    if (cg->side == 0){
        LocalParticle_kill_particle(part, XC_ERR_NOT_IMPLEMENTED);
    }

    SRotation_single_particle(part, cg->sin_z, cg->cos_z);
    part_x = LocalParticle_get_x(part);
#ifdef XCOLL_USE_EXACT
    part_tan_x = LocalParticle_get_exact_xp(part);
#else
    part_tan_x = LocalParticle_get_xp(part);
#endif
    part_y = LocalParticle_get_y(part);
#ifdef XCOLL_USE_EXACT
    part_tan_y = LocalParticle_get_exact_yp(part);
#else
    part_tan_y = LocalParticle_get_yp(part);
#endif
    s = get_s_of_first_crossing_with_vlimit(part_x, part_tan_x, part_y, part_tan_y, cg->segments, 4, -cg->height/2, cg->height/2);
    if (s < S_MAX){
        is_hit = cg->side;
    } else {
        // No hit, rotate back to lab frame
        SRotation_single_particle(part, -cg->sin_z, cg->cos_z);
    }

    // Drift to the impact position or end, and move to jaw frame if relevant
    if (is_hit != 0){
        // Move to the impact position
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, s);
#else
        Drift_single_particle_expanded(part, s);
#endif
        // Shift the reference frame to the jaw corner LU
        XYShift_single_particle(part, cg->jaw_U, 0);
        LocalParticle_add_to_s(part, -cg->length/2*(1 - cg->cos_y));
        // Rotate the reference frame to tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), asin(cg->sin_y));
        LocalParticle_set_s(part, new_s);
        if (cg->side == 1){
            if (cg->record_touches){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_L);
            }
        } else {
            // Mirror x
            LocalParticle_scale_x(part, -1);
            LocalParticle_scale_px(part, -1);
            cg->bending_radius = -cg->bending_radius;
            cg->bending_angle = -cg->bending_angle;
            if (cg->record_touches){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_R);
            }
        }

    } else {
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, cg->length);
#else
        Drift_single_particle_expanded(part, cg->length);
#endif
    }

    return is_hit;
}


/*gpufun*/
void hit_crystal_transform_back(int8_t is_hit, LocalParticle* part, CrystalGeometry cg){
    if (is_hit != 0){
        if (cg->side == -1){
            // Mirror back
            LocalParticle_scale_x(part, -1);
            LocalParticle_scale_px(part, -1);
            cg->bending_radius = -cg->bending_radius;
            cg->bending_angle = -cg->bending_angle;
        }
        // Rotate back from tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), -asin(cg->sin_y));
        LocalParticle_set_s(part, new_s);
        // Shift the reference frame back from jaw corner U
        XYShift_single_particle(part, -cg->jaw_U, 0);
        LocalParticle_add_to_s(part, cg->length/2*(1 - cg->cos_y));
        // If particle survived, drift to end of element
        if (LocalParticle_get_state(part) > 0){
            if (cg->record_touches){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_EXIT_JAW);
            }
#ifdef XCOLL_USE_EXACT
            Drift_single_particle_exact(part, cg->length - LocalParticle_get_s(part));
#else
            Drift_single_particle_expanded(part, cg->length - LocalParticle_get_s(part));
#endif
        }
        SRotation_single_particle(part, -cg->sin_z, cg->cos_z);
    }
}


#endif /* XCOLL_CRY_GEOM_H */
