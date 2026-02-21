// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_CRY_GEOM_H
#define XCOLL_CRY_GEOM_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU


#include <xtrack/headers/track.h>
#include <xcoll/headers/particle_states.h>
#include <xcoll/scattering_routines/geometry/get_s.h>
#include <xcoll/scattering_routines/geometry/rotation.h>


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
    double bending_angle;
    double width;
    double height;
    double miscut_angle;
    double s_B;    // Bend centre
    double x_B;
    double s_P;    // Miscut centre
    double x_P;
    double t_VImax;
    // Segments
    Segment* segments;
    // Impact table
    InteractionRecordData record;
    RecordIndex record_index;
    int8_t record_impacts;
    int8_t record_exits;
} CrystalGeometry_;
typedef CrystalGeometry_* CrystalGeometry;


// This function checks if a particle hits a jaw (and which).
// Return value: 0 (no hit), 1 (hit on left jaw), -1 (hit on right jaw).
/*gpufun*/
int8_t hit_crystal_check(LocalParticle* part, CrystalGeometry restrict cg){
    double part_x, part_tan_x, part_y, part_tan_y;
    double s = 1.e21;

    // Crystal should be single-sided
    if (cg->side == 0){
        LocalParticle_kill_particle(part, XC_ERR_NOT_IMPLEMENTED);
        return 0;
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

    // Rotate back to lab frame
    SRotation_single_particle(part, -cg->sin_z, cg->cos_z);

    if (s < S_MAX){
        return cg->side;

    } else {
        return 0;
    }
}


// This function checks if a particle hits a jaw (and which).
// Return value: 0 (no hit), 1 (hit on left jaw), -1 (hit on right jaw).
// Furthermore, the particle is moved to the location where it hits the jaw (drifted to the end if no hit),
//              and transformed to the reference frame of that jaw.
/*gpufun*/
int8_t hit_crystal_check_and_transform(LocalParticle* part, CrystalGeometry restrict cg){
    double part_x, part_tan_x, part_y, part_tan_y;
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
        // Hit: Drift to the impact position, and move to jaw frame if relevant
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, s);
#else
        Drift_single_particle_expanded(part, s);
#endif
        // Shift the reference frame to the upstream jaw corner (for a crystal, this is always at s=0)
        XYShift_single_particle(part, cg->jaw_U, 0);
        // Rotate the reference frame to tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), asin(cg->sin_y));
        LocalParticle_set_s(part, new_s);
        if (cg->side == 1){
            if (cg->record_impacts){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_L);
            }

        } else {
            // Mirror x
            LocalParticle_scale_x(part, -1);
#ifdef XCOLL_USE_EXACT
            LocalParticle_scale_exact_xp(part, -1);
#else
            LocalParticle_scale_xp(part, -1);
#endif
            if (cg->record_impacts){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_R);
            }
        }
        return cg->side;

    } else {
        // No hit, rotate back to lab frame
        SRotation_single_particle(part, -cg->sin_z, cg->cos_z);
        // Drift to end
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, cg->length);
#else
        Drift_single_particle_expanded(part, cg->length);
#endif
        return 0;
    }
}


// Return to start position after having logged the impact.
/*gpufun*/
void hit_crystal_return(int8_t is_hit, LocalParticle* part, CrystalGeometry restrict cg){
    if (is_hit != 0){
        if (cg->side == -1){
            // Mirror back
            LocalParticle_scale_x(part, -1);
#ifdef XCOLL_USE_EXACT
            LocalParticle_scale_exact_xp(part, -1);
#else
            LocalParticle_scale_xp(part, -1);
#endif
        }
        // Rotate back from tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), -asin(cg->sin_y));
        LocalParticle_set_s(part, new_s);
        // Shift the reference frame back from jaw corner U
        XYShift_single_particle(part, -cg->jaw_U, 0);
        // Drift back to start of element
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, -LocalParticle_get_s(part));
#else
        Drift_single_particle_expanded(part, -LocalParticle_get_s(part));
#endif
        SRotation_single_particle(part, -cg->sin_z, cg->cos_z);
    }
}


/*gpufun*/
void hit_crystal_transform_back(int8_t is_hit, LocalParticle* part, CrystalGeometry restrict cg){
    if (is_hit != 0){
        if (LocalParticle_get_state(part) > 0){
            if (cg->record_exits){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_EXIT_JAW);
            }
        }
        if (cg->side == -1){
            // Mirror back
            LocalParticle_scale_x(part, -1);
#ifdef XCOLL_USE_EXACT
            LocalParticle_scale_exact_xp(part, -1);
#else
            LocalParticle_scale_xp(part, -1);
#endif
        }
        // Rotate back from tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), -asin(cg->sin_y));
        LocalParticle_set_s(part, new_s);
        // Shift the reference frame back from jaw corner U
        XYShift_single_particle(part, -cg->jaw_U, 0);
        // If particle survived, drift to end of element
        if (LocalParticle_get_state(part) > 0){
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
