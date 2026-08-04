// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_FLUKA_CRYSTAL_H
#define XCOLL_FLUKA_CRYSTAL_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#endif /* XO_CONTEXT_CPU */

#include "xobjects/headers/common.h"
#include "xtrack/headers/checks.h"
#include "xtrack/headers/track.h"
#include "xtrack/beam_elements/elements_src/track_drift.h"
#include "xcoll/headers/checks.h"
#include "xcoll/scattering_routines/geometry/objects.h"
#include "xcoll/scattering_routines/geometry/crystal_geometry.h"


GPUFUN
int8_t FlukaCrystalData_get_record_impacts(FlukaCrystalData el) {
    return FlukaCrystalData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t FlukaCrystalData_get_record_exits(FlukaCrystalData el) {
    return (FlukaCrystalData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t FlukaCrystalData_get_record_scatterings(FlukaCrystalData el) {
    return (FlukaCrystalData_get__record_interactions(el) >> 2) % 2;
}


GPUFUN
void FlukaCrystal_init_geometry(FlukaCrystalData el, LocalParticle* part0,
                                CrystalGeometry RESTRICT cg) {
    // Dimensions
    cg->length = FlukaCrystalData_get_length(el);
    cg->side   = FlukaCrystalData_get__side(el);
    cg->bending_radius = FlukaCrystalData_get__bending_radius(el);
    cg->bending_angle  = FlukaCrystalData_get__bending_angle(el);
    cg->width  = FlukaCrystalData_get__width(el);
    cg->height = FlukaCrystalData_get__height(el);
    if (cg->width<=0) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Crystal width should be positive!");
        return;
    }
    if (cg->height<=0) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Crystal height should be positive!");
        return;
    }
    if (cg->bending_radius==0) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Crystal bending radius should be non-zero!");
        return;
    }
    if (cg->side != 1 && cg->side != -1) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Crystal side should be either 1 or -1!");
        return;
    }
    // Jaw corners (with tilts)
    cg->jaw_U = FlukaCrystalData_get__jaw_U(el);
    // Jaw angles
    cg->sin_z = FlukaCrystalData_get__sin_z(el);
    cg->cos_z = FlukaCrystalData_get__cos_z(el);
    // Jaw tilts
    cg->sin_y = FlukaCrystalData_get__sin_y(el);
    cg->cos_y = FlukaCrystalData_get__cos_y(el);
    // Segments
    int8_t status = 0;
    double jaw;
    if (cg->side == 1) {
        jaw = cg->jaw_U;
    } else if (cg->side == -1) {
        jaw = cg->jaw_U - cg->width; // To ensure that jaw_U is the inner corner
    }
    status = create_crystal(cg->segments, cg->bending_radius, cg->width,
                            cg->length, jaw, cg->sin_y, cg->cos_y);
    if (status != 0) {
        xcoll_kill_with_message(part0, XC_ERR_FAILED_GEOMETRY,
                                "Failed creating crystal geometry!");
        return;
    }
    // Impact table
    cg->record = FlukaCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record) {
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = FlukaCrystalData_get_record_impacts(el);
        cg->record_exits = FlukaCrystalData_get_record_exits(el);
    }
    // Not needed, set to zero
    cg->miscut_angle = 0;
    cg->s_P = 0;
    cg->x_P = 0;
    cg->t_VImax = 0;
}


GPUFUN
void FlukaCrystal_track_local_particle(FlukaCrystalData el, LocalParticle* part0) {
    int8_t active = FlukaCrystalData_get_active(el);
    active       *= FlukaCrystalData_get__tracking(el);
    double const length = FlukaCrystalData_get_length(el);
    double const length_front = FlukaCrystalData_get_length_front(el);
    double const length_back = FlukaCrystalData_get_length_back(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    CrystalGeometry_ cg_;
    CrystalGeometry cg = &cg_; // pointer
    if (active) {
        FlukaCrystal_init_geometry(el, part0, cg);
    }

    START_PER_PARTICLE_BLOCK(part0, part);
        if (!active) {
            // Drift full length
            Drift_single_particle(part, length+length_front+length_back);

        } else {
            // Check collimator initialisation
            int8_t is_valid = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_valid) {
                // Move to start of crystal
#ifdef XCOLL_USE_EXACT
                Drift_single_particle_exact(part, length_front);
#else
                Drift_single_particle_expanded(part, length_front);
#endif

                // Store s-location of start of crystal
                double const s_coll = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Check if hit on jaws
                int8_t is_hit;
                if (cg->record_impacts) {
                    // Check with transformation (to log impacts).
                    // Particle will be at impact position or at end.
                    is_hit = hit_crystal_check_and_transform(part, cg);
                    if (is_hit != 0) {
                        // Particle needs to return to start position.
                        hit_crystal_return(is_hit, part, cg);
                    }

                } else {
                    // Only check. Particle will be at start position.
                    is_hit = hit_crystal_check(part, cg);
                    if (is_hit == 0) {
                        // Drift to end.
#ifdef XCOLL_USE_EXACT
                        Drift_single_particle_exact(part, cg->length);
#else
                        Drift_single_particle_expanded(part, cg->length);
#endif
                    }
                }

                if (is_hit == 0) {
                    // Drift to end.
                    double length_back = FlukaCrystalData_get_length_back(el);
#ifdef XCOLL_USE_EXACT
                    Drift_single_particle_exact(part, length_back);
#else
                    Drift_single_particle_expanded(part, length_back);
#endif

                } else {
                    // Mark for FLUKA processing.
                    if (LocalParticle_get_state(part) == XC_SECONDARY_PARTICLE) {
                        LocalParticle_set_state(part, XC_HIT_ON_FLUKA_SEC);
                    } else {
                        LocalParticle_set_state(part, XC_HIT_ON_FLUKA);
                    }
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
    END_PER_PARTICLE_BLOCK;
}

#endif /* XCOLL_FLUKA_CRYSTAL_H */
