// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_TRANSPARENT_CRYSTAL_H
#define XCOLL_TRANSPARENT_CRYSTAL_H

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
int8_t TransparentCrystalData_get_record_impacts(TransparentCrystalData el) {
    return TransparentCrystalData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t TransparentCrystalData_get_record_exits(TransparentCrystalData el) {
    return (TransparentCrystalData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t TransparentCrystalData_get_record_scatterings(TransparentCrystalData el) {
    return (TransparentCrystalData_get__record_interactions(el) >> 2) % 2;
}


GPUFUN
void TransparentCrystal_init_geometry(TransparentCrystalData el, LocalParticle* part0,
                                      CrystalGeometry RESTRICT cg) {
    // Dimensions
    cg->length = TransparentCrystalData_get_length(el);
    cg->side   = TransparentCrystalData_get__side(el);
    cg->bending_radius = TransparentCrystalData_get__bending_radius(el);
    cg->bending_angle  = TransparentCrystalData_get__bending_angle(el);
    cg->width  = TransparentCrystalData_get__width(el);
    cg->height = TransparentCrystalData_get__height(el);
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
    cg->jaw_U = TransparentCrystalData_get__jaw_U(el);
    // Jaw angles
    cg->sin_z = TransparentCrystalData_get__sin_z(el);
    cg->cos_z = TransparentCrystalData_get__cos_z(el);
    // Jaw tilts
    cg->sin_y = TransparentCrystalData_get__sin_y(el);
    cg->cos_y = TransparentCrystalData_get__cos_y(el);
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
    cg->record = TransparentCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record) {
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = TransparentCrystalData_get_record_impacts(el);
        cg->record_exits = TransparentCrystalData_get_record_exits(el);
    }
    // Not needed, set to zero
    cg->miscut_angle = 0;
    cg->s_P = 0;
    cg->x_P = 0;
    cg->t_VImax = 0;
}


GPUFUN
void TransparentCrystal_track_local_particle(TransparentCrystalData el, LocalParticle* part0) {
    int8_t active = TransparentCrystalData_get_active(el);
    active       *= TransparentCrystalData_get__tracking(el);
    double const length = TransparentCrystalData_get_length(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    CrystalGeometry_ cg_;
    CrystalGeometry cg = &cg_; // pointer
    if (active) {
        TransparentCrystal_init_geometry(el, part0, cg);
    }

    START_PER_PARTICLE_BLOCK(part0, part);
        if (!active) {
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_valid = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_valid) {
                // Store s-location of start of crystal
                double const s_coll = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Check if hit on jaws
                int8_t is_hit = hit_crystal_check_and_transform(part, cg);

                // Transform back to the lab frame
                hit_crystal_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);
            }
        }
    END_PER_PARTICLE_BLOCK;
}

#endif /* XCOLL_TRANSPARENT_CRYSTAL_H */
