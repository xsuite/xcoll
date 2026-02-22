// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEANT4_TIP_H
#define XCOLL_GEANT4_TIP_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU

#include <xobjects/headers/common.h>
#include <xtrack/headers/track.h>
#include <xtrack/headers/checks.h>
#include <xtrack/beam_elements/elements_src/track_drift.h>
#include <xcoll/headers/checks.h>
#include <xcoll/lib/particle_states.h>      // auto-generated from xcoll/headers/particle_states.py
#include <xcoll/scattering_routines/geometry/objects.h>
#include <xcoll/scattering_routines/geometry/collimator_geometry.h>


GPUFUN
int8_t Geant4CollimatorTipData_get_record_impacts(Geant4CollimatorTipData el){
    return Geant4CollimatorTipData_get__record_interactions(el) % 2;
}

GPUFUN
int8_t Geant4CollimatorTipData_get_record_exits(Geant4CollimatorTipData el){
    return (Geant4CollimatorTipData_get__record_interactions(el) >> 1) % 2;
}

GPUFUN
int8_t Geant4CollimatorTipData_get_record_scatterings(Geant4CollimatorTipData el){
    return (Geant4CollimatorTipData_get__record_interactions(el) >> 2) % 2;
}


GPUFUN
CollimatorGeometry Geant4CollimatorTip_init_geometry(Geant4CollimatorTipData el, LocalParticle* part0, int8_t active){
    CollimatorGeometry cg = (CollimatorGeometry) malloc(sizeof(CollimatorGeometry_));
    if (active){ // This is needed in order to avoid that the initialisation is called during a twiss!
        // Jaw corners (with tilts)
        cg->jaw_LU = Geant4CollimatorTipData_get__jaw_LU(el);
        cg->jaw_RU = Geant4CollimatorTipData_get__jaw_RU(el);
        // Get angles of jaws
        cg->sin_zL = Geant4CollimatorTipData_get__sin_zL(el);
        cg->cos_zL = Geant4CollimatorTipData_get__cos_zL(el);
        cg->sin_zR = Geant4CollimatorTipData_get__sin_zR(el);
        cg->cos_zR = Geant4CollimatorTipData_get__cos_zR(el);
        cg->sin_zDiff = Geant4CollimatorTipData_get__sin_zDiff(el);
        cg->cos_zDiff = Geant4CollimatorTipData_get__cos_zDiff(el);
        cg->jaws_parallel = Geant4CollimatorTipData_get__jaws_parallel(el);
        // Tilts
        cg->sin_yL = Geant4CollimatorTipData_get__sin_yL(el);
        cg->cos_yL = Geant4CollimatorTipData_get__cos_yL(el);
        cg->sin_yR = Geant4CollimatorTipData_get__sin_yR(el);
        cg->cos_yR = Geant4CollimatorTipData_get__cos_yR(el);
        // Length and segments
        cg->length = Geant4CollimatorTipData_get_length(el);
        cg->side   = Geant4CollimatorTipData_get__side(el);
        double s_U, s_D, x_D;
        if (cg->side != -1){
            s_U = cg->length/2 * (1-cg->cos_yL);
            s_D = cg->length/2 * (1+cg->cos_yL);
            x_D = Geant4CollimatorTipData_get__jaw_LD(el);
            cg->segments_L = create_jaw(s_U, cg->jaw_LU, s_D, x_D, cg->sin_yL/cg->cos_yL, 1);
        }
        if (cg->side != 1){
            s_U = cg->length/2 * (1-cg->cos_yR);
            s_D = cg->length/2 * (1+cg->cos_yR);
            x_D = Geant4CollimatorTipData_get__jaw_RD(el);
            cg->segments_R = create_jaw(s_U, cg->jaw_RU, s_D, x_D, cg->sin_yR/cg->cos_yR, -1);
        }
        // Impact table
        cg->record = Geant4CollimatorTipData_getp_internal_record(el, part0);
        cg->record_index = NULL;
        cg->record_impacts = 0;
        cg->record_exits = 0;
        if (cg->record){
            cg->record_index = InteractionRecordData_getp__index(cg->record);
            cg->record_impacts = Geant4CollimatorTipData_get_record_impacts(el);
            cg->record_exits = Geant4CollimatorTipData_get_record_exits(el);
        }
    }

    return cg;
}

GPUFUN
void Geant4CollimatorTip_free(CollimatorGeometry restrict cg, int8_t active){
    if (active){
        if (cg->side != -1){
            destroy_jaw(cg->segments_L);
        }
        if (cg->side != 1){
            destroy_jaw(cg->segments_R);
        }
    }
    free(cg);
}


GPUFUN
void Geant4CollimatorTip_track_local_particle(Geant4CollimatorTipData el, LocalParticle* part0){

    // Collimator active and length
    int8_t active = Geant4CollimatorTipData_get_active(el);
    active       *= Geant4CollimatorTipData_get__tracking(el);
    double const length = Geant4CollimatorTipData_get_length(el);

    // Get geometry
    CollimatorGeometry cg = Geant4CollimatorTip_init_geometry(el, part0, active);

    START_PER_PARTICLE_BLOCK(part0, part);
        if (!active){
            // Drift full length (use global setting for expanded vs exact drift)
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_tracking) {
                LocalParticle_set_s(part, 0);

                // Check if hit on jaws
                int8_t is_hit = hit_jaws_check_and_transform(part, cg);
                (void) is_hit;
            }
        }
    END_PER_PARTICLE_BLOCK;
    Geant4CollimatorTip_free(cg, active);
}

#endif /* XCOLL_GEANT4_TIP_H */
