// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_FLUKA_COLLIMATOR_H
#define XCOLL_FLUKA_COLLIMATOR_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#endif /* XO_CONTEXT_CPU */

#include "xobjects/headers/common.h"
#include "xtrack/headers/checks.h"
#include "xtrack/headers/track.h"
#include "xtrack/beam_elements/elements_src/track_drift.h"
#include "xcoll/headers/checks.h"
#include "xcoll/scattering_routines/geometry/objects.h"
#include "xcoll/scattering_routines/geometry/collimator_geometry.h"


GPUFUN
int8_t FlukaCollimatorData_get_record_impacts(FlukaCollimatorData el) {
    return FlukaCollimatorData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t FlukaCollimatorData_get_record_exits(FlukaCollimatorData el) {
    return (FlukaCollimatorData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t FlukaCollimatorData_get_record_scatterings(FlukaCollimatorData el) {
    return (FlukaCollimatorData_get__record_interactions(el) >> 2) % 2;
}


GPUFUN
void FlukaCollimator_init_geometry(FlukaCollimatorData el, LocalParticle* part0,
                                   CollimatorGeometry RESTRICT cg) {
    // Dimensions
    cg->length  = FlukaCollimatorData_get_length(el);
    cg->length += FlukaCollimatorData_get_length_front(el);
    cg->length += FlukaCollimatorData_get_length_back(el);
    cg->side   = FlukaCollimatorData_get__side(el);
    if (cg->side != 1 && cg->side != -1 && cg->side != 0) {
        xcoll_kill_with_message(part0, XC_ERR_INVALID_XOFIELD,
                                "Collimator side should be either 1, -1, or 0!");
        return;
    }
    // Jaw corners (with tilts)
    cg->jaw_LU = FlukaCollimatorData_get__jaw_LU(el);
    cg->jaw_RU = FlukaCollimatorData_get__jaw_RU(el);
    // Jaw angles
    cg->sin_zL = FlukaCollimatorData_get__sin_zL(el);
    cg->cos_zL = FlukaCollimatorData_get__cos_zL(el);
    cg->sin_zR = FlukaCollimatorData_get__sin_zR(el);
    cg->cos_zR = FlukaCollimatorData_get__cos_zR(el);
    cg->sin_zDiff = FlukaCollimatorData_get__sin_zDiff(el);
    cg->cos_zDiff = FlukaCollimatorData_get__cos_zDiff(el);
    cg->jaws_parallel = FlukaCollimatorData_get__jaws_parallel(el);
    // Jaw tilts
    cg->sin_yL = FlukaCollimatorData_get__sin_yL(el);
    cg->cos_yL = FlukaCollimatorData_get__cos_yL(el);
    cg->sin_yR = FlukaCollimatorData_get__sin_yR(el);
    cg->cos_yR = FlukaCollimatorData_get__cos_yR(el);
    // Segments
    int8_t status = 0;
    double s_U, s_D, x_D;
    if (cg->side != -1) {
        s_U = cg->length/2 * (1-cg->cos_yL);
        s_D = cg->length/2 * (1+cg->cos_yL);
        x_D = FlukaCollimatorData_get__jaw_LD(el);
        status = create_jaw(cg->segments_L, s_U, cg->jaw_LU, s_D, x_D, cg->sin_yL/cg->cos_yL, 1);
    }
    if (cg->side != 1) {
        s_U = cg->length/2 * (1-cg->cos_yR);
        s_D = cg->length/2 * (1+cg->cos_yR);
        x_D = FlukaCollimatorData_get__jaw_RD(el);
        status = create_jaw(cg->segments_R, s_U, cg->jaw_RU, s_D, x_D, cg->sin_yR/cg->cos_yR, -1);
    }
    if (status != 0) {
        xcoll_kill_with_message(part0, XC_ERR_FAILED_GEOMETRY,
                                "Failed creating collimator geometry!");
        return;
    }
    // Impact table
    cg->record = FlukaCollimatorData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record) {
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = FlukaCollimatorData_get_record_impacts(el);
        cg->record_exits = FlukaCollimatorData_get_record_exits(el);
    }
}


GPUFUN
void FlukaCollimator_track_local_particle(FlukaCollimatorData el, LocalParticle* part0) {
    int8_t active = FlukaCollimatorData_get_active(el);
    active       *= FlukaCollimatorData_get__tracking(el);
    double const length = FlukaCollimatorData_get_length(el);
    double const length_front = FlukaCollimatorData_get_length_front(el);
    double const length_back = FlukaCollimatorData_length_back(el);

    // Initialise collimator data
    // TODO: we want this to happen before tracking (instead of every turn), as a separate kernel
    CollimatorGeometry_ cg_;
    CollimatorGeometry cg = &cg_; // pointer
    if (active) {
        FlukaCollimator_init_geometry(el, part0, cg);
    }

    START_PER_PARTICLE_BLOCK(part0, part);
        if (!active) {
            // Drift full length
            Drift_single_particle(part, length+length_front+length_back);

        } else {
            // Check collimator initialisation
            int8_t is_valid = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_valid) {
                // Store s-location of start of collimator
                double const s_coll = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Check if hit on jaws
                int8_t is_hit;
                if (cg->record_impacts) {
                    // Check with transformation (to log impacts).
                    // Particle will be at impact position or at end.
                    is_hit = hit_jaws_check_and_transform(part, cg);
                    if (is_hit != 0) {
                        // Particle needs to return to start position.
                        hit_jaws_return(is_hit, part, cg);
                    }

                } else {
                    // Only check. Particle will be at start position.
                    is_hit = hit_jaws_check(part, cg);
                    if (is_hit == 0) {
                        // Drift to end.
#ifdef XCOLL_USE_EXACT
                        Drift_single_particle_exact(part, cg->length);
#else
                        Drift_single_particle_expanded(part, cg->length);
#endif
                    }
                }

                if (is_hit != 0) {
                    // Mark for FLUKA processing.
                    if (LocalParticle_get_state(part) == XC_SECONDARY_PARTICLE) {
                        LocalParticle_set_state(part, XC_HIT_ON_FLUKA_SEC);
                    } else {
                        LocalParticle_set_state(part, XC_HIT_ON_FLUKA);
                    }
                }
                LocalParticle_add_to_s(part, s_coll);
            }
        }
    END_PER_PARTICLE_BLOCK;
}

#endif /* XCOLL_FLUKA_COLLIMATOR_H */
