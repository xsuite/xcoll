// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_ABSORBER_H
#define XCOLL_ABSORBER_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU

#include "xobjects/headers/common.h"
#include "xtrack/headers/track.h"
#include "xtrack/headers/checks.h"
#include "xtrack/beam_elements/elements_src/track_drift.h"
#include "xcoll/generated_src/particle_states.h"      // from xcoll/headers/particle_states.py
#include "xcoll/generated_src/interaction_types.h"    // from xcoll/interaction_record/interaction_types.py
#include "xcoll/interaction_record/interaction_record_src/interaction_record.h"
#include "xcoll/scattering_routines/geometry/objects.h"
#include "xcoll/scattering_routines/geometry/collimator_geometry.h"


GPUFUN
int8_t BlackAbsorberData_get_record_impacts(BlackAbsorberData el){
    return BlackAbsorberData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t BlackAbsorberData_get_record_exits(BlackAbsorberData el){
    return (BlackAbsorberData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t BlackAbsorberData_get_record_scatterings(BlackAbsorberData el){
    return (BlackAbsorberData_get__record_interactions(el) >> 2) % 2;
}


GPUFUN
CollimatorGeometry BlackAbsorber_init_geometry(BlackAbsorberData el, LocalParticle* part0){
    CollimatorGeometry cg = (CollimatorGeometry) malloc(sizeof(CollimatorGeometry_));
    // Jaw corners (with tilts)
    cg->jaw_LU = BlackAbsorberData_get__jaw_LU(el);
    cg->jaw_RU = BlackAbsorberData_get__jaw_RU(el);
    // Get angles of jaws
    cg->sin_zL = BlackAbsorberData_get__sin_zL(el);
    cg->cos_zL = BlackAbsorberData_get__cos_zL(el);
    cg->sin_zR = BlackAbsorberData_get__sin_zR(el);
    cg->cos_zR = BlackAbsorberData_get__cos_zR(el);
    cg->sin_zDiff = BlackAbsorberData_get__sin_zDiff(el);
    cg->cos_zDiff = BlackAbsorberData_get__cos_zDiff(el);
    cg->jaws_parallel = BlackAbsorberData_get__jaws_parallel(el);
    // Tilts
    cg->sin_yL = BlackAbsorberData_get__sin_yL(el);
    cg->cos_yL = BlackAbsorberData_get__cos_yL(el);
    cg->sin_yR = BlackAbsorberData_get__sin_yR(el);
    cg->cos_yR = BlackAbsorberData_get__cos_yR(el);
    // Length and segments
    cg->length = BlackAbsorberData_get_length(el);
    cg->side   = BlackAbsorberData_get__side(el);
    double s_U, s_D, x_D;
    if (cg->side != -1){
        s_U = cg->length/2 * (1-cg->cos_yL);
        s_D = cg->length/2 * (1+cg->cos_yL);
        x_D = BlackAbsorberData_get__jaw_LD(el);
        cg->segments_L = create_jaw(s_U, cg->jaw_LU, s_D, x_D, cg->sin_yL/cg->cos_yL, 1);
    }
    if (cg->side != 1){
        s_U = cg->length/2 * (1-cg->cos_yR);
        s_D = cg->length/2 * (1+cg->cos_yR);
        x_D = BlackAbsorberData_get__jaw_RD(el);
        cg->segments_R = create_jaw(s_U, cg->jaw_RU, s_D, x_D, cg->sin_yR/cg->cos_yR, -1);
    }
    // Impact table
    cg->record = BlackAbsorberData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record){
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = BlackAbsorberData_get_record_impacts(el);
        cg->record_exits = BlackAbsorberData_get_record_exits(el);
    }
    return cg;
}

GPUFUN
void BlackAbsorber_free(CollimatorGeometry restrict cg){
    if (cg->side != -1){
        destroy_jaw(cg->segments_L);
    }
    if (cg->side != 1){
        destroy_jaw(cg->segments_R);
    }
    free(cg);
}


GPUFUN
void BlackAbsorber_track_local_particle(BlackAbsorberData el, LocalParticle* part0){
    int8_t active = BlackAbsorberData_get_active(el);
    active       *= BlackAbsorberData_get__tracking(el);
    double const length = BlackAbsorberData_get_length(el);

    // Get geometry
    CollimatorGeometry cg;
    int8_t record_scatterings;
    if (active){
        cg = BlackAbsorber_init_geometry(el, part0);
        record_scatterings = BlackAbsorberData_get_record_scatterings(el);
    }

    START_PER_PARTICLE_BLOCK(part0, part);
        if (!active){
            // Drift full length (use global setting for expanded vs exact drift)
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_tracking) {
                // Store s-location of start of collimator
                double s_coll = LocalParticle_get_s(part);
                LocalParticle_set_s(part, 0);

                // Check if hit on jaws
                int8_t is_hit = hit_jaws_check_and_transform(part, cg);

                if (is_hit != 0){
                    LocalParticle_set_state(part, XC_LOST_ON_BLACK_ABSORBER);
                    if (record_scatterings) {
                        InteractionRecordData_log(cg->record, cg->record_index, part, XC_ABSORBED);  // In coll jaw reference frame
                    }
                }

                // Transform back to the lab frame
                hit_jaws_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);
            }
        }
    END_PER_PARTICLE_BLOCK;
    if (active){
        BlackAbsorber_free(cg);
    }
}

#endif /* XCOLL_ABSORBER_H */
