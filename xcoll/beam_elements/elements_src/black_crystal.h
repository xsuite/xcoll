// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_BLACK_CRYSTAL_H
#define XCOLL_BLACK_CRYSTAL_H

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
#include "xcoll/interaction_record/interaction_record_src/interaction_record.h"


GPUFUN
int8_t BlackCrystalData_get_record_impacts(BlackCrystalData el) {
    return BlackCrystalData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t BlackCrystalData_get_record_exits(BlackCrystalData el) {
    return (BlackCrystalData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t BlackCrystalData_get_record_scatterings(BlackCrystalData el) {
    return (BlackCrystalData_get__record_interactions(el) >> 2) % 2;
}


GPUFUN
void BlackCrystal_init_geometry(BlackCrystalData el, LocalParticle* part0,
                                CrystalGeometry RESTRICT cg) {
    // Dimensions
    cg->length = BlackCrystalData_get_length(el);
    cg->side   = BlackCrystalData_get__side(el);
    cg->bending_radius = BlackCrystalData_get__bending_radius(el);
    cg->bending_angle  = BlackCrystalData_get__bending_angle(el);
    cg->width  = BlackCrystalData_get__width(el);
    cg->height = BlackCrystalData_get__height(el);
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
    cg->jaw_U = BlackCrystalData_get__jaw_U(el);
    // Jaw angles
    cg->sin_z = BlackCrystalData_get__sin_z(el);
    cg->cos_z = BlackCrystalData_get__cos_z(el);
    // Jaw tilts
    cg->sin_y = BlackCrystalData_get__sin_y(el);
    cg->cos_y = BlackCrystalData_get__cos_y(el);
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
    cg->record = BlackCrystalData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record) {
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = BlackCrystalData_get_record_impacts(el);
        cg->record_exits = BlackCrystalData_get_record_exits(el);
    }
    // Not needed, set to zero
    cg->miscut_angle = 0;
    cg->s_P = 0;
    cg->x_P = 0;
    cg->t_VImax = 0;
}


GPUFUN
void BlackCrystal_track_local_particle(BlackCrystalData el, LocalParticle* part0) {
    int8_t active = BlackCrystalData_get_active(el);
    active       *= BlackCrystalData_get__tracking(el);
    double const length = BlackCrystalData_get_length(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    CrystalGeometry_ cg_;
    CrystalGeometry cg = &cg_; // pointer
    int8_t record_scatterings;
    if (active) {
        BlackCrystal_init_geometry(el, part0, cg);
        record_scatterings = BlackCrystalData_get_record_scatterings(el);
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

                if (is_hit != 0) {
                    if (LocalParticle_get_state(part) == XC_SECONDARY_PARTICLE) {
                        LocalParticle_set_state(part, XC_LOST_ON_MATERIAL_SEC);
                    } else {
                        LocalParticle_set_state(part, XC_LOST_ON_MATERIAL);
                    }
                    if (record_scatterings) {
                        // In jaw reference frame
                        InteractionRecordData_log(cg->record, cg->record_index,
                                                  part, XC_ABSORBED, is_hit);
                    }
                }

                // Transform back to the lab frame
                hit_crystal_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);
            }
        }
    END_PER_PARTICLE_BLOCK;
}

#endif /* XCOLL_BLACK_CRYSTAL_H */
