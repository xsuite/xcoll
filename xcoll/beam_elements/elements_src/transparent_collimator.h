// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_TRANSPARENT_COLL_H
#define XCOLL_TRANSPARENT_COLL_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU


/*gpufun*/
int8_t TransparentCollimatorData_get_record_impacts(TransparentCollimatorData el){
    return TransparentCollimatorData_get__record_interactions(el) % 2;
}

/*gpufun*/
int8_t TransparentCollimatorData_get_record_exits(TransparentCollimatorData el){
    return (TransparentCollimatorData_get__record_interactions(el) >> 1) % 2;
}

/*gpufun*/
int8_t TransparentCollimatorData_get_record_scatterings(TransparentCollimatorData el){
    return (TransparentCollimatorData_get__record_interactions(el) >> 2) % 2;
}


/*gpufun*/
CollimatorGeometry TransparentCollimator_init_geometry(TransparentCollimatorData el, LocalParticle* part0){
    CollimatorGeometry cg = (CollimatorGeometry) malloc(sizeof(CollimatorGeometry_));
    // Jaw corners (with tilts)
    cg->jaw_LU = TransparentCollimatorData_get__jaw_LU(el);
    cg->jaw_RU = TransparentCollimatorData_get__jaw_RU(el);
    // Get angles of jaws
    cg->sin_zL = TransparentCollimatorData_get__sin_zL(el);
    cg->cos_zL = TransparentCollimatorData_get__cos_zL(el);
    cg->sin_zR = TransparentCollimatorData_get__sin_zR(el);
    cg->cos_zR = TransparentCollimatorData_get__cos_zR(el);
    cg->sin_zDiff = TransparentCollimatorData_get__sin_zDiff(el);
    cg->cos_zDiff = TransparentCollimatorData_get__cos_zDiff(el);
    cg->jaws_parallel = TransparentCollimatorData_get__jaws_parallel(el);
    // Tilts
    cg->sin_yL = TransparentCollimatorData_get__sin_yL(el);
    cg->cos_yL = TransparentCollimatorData_get__cos_yL(el);
    cg->sin_yR = TransparentCollimatorData_get__sin_yR(el);
    cg->cos_yR = TransparentCollimatorData_get__cos_yR(el);
    // Length and segments
    cg->length = TransparentCollimatorData_get_length(el);
    cg->side   = TransparentCollimatorData_get__side(el);
    double s_U, s_D, x_D;
    if (cg->side != -1){
        s_U = cg->length/2 * (1-cg->cos_yL);
        s_D = cg->length/2 * (1+cg->cos_yL);
        x_D = TransparentCollimatorData_get__jaw_LD(el);
        cg->segments_L = create_jaw(s_U, cg->jaw_LU, s_D, x_D, cg->sin_yL/cg->cos_yL, 1);
    }
    if (cg->side != 1){
        s_U = cg->length/2 * (1-cg->cos_yR);
        s_D = cg->length/2 * (1+cg->cos_yR);
        x_D = TransparentCollimatorData_get__jaw_RD(el);
        cg->segments_R = create_jaw(s_U, cg->jaw_RU, s_D, x_D, cg->sin_yR/cg->cos_yR, -1);
    }
    // Impact table
    cg->record = TransparentCollimatorData_getp_internal_record(el, part0);
    cg->record_index = NULL;
    cg->record_impacts = 0;
    cg->record_exits = 0;
    if (cg->record){
        cg->record_index = InteractionRecordData_getp__index(cg->record);
        cg->record_impacts = TransparentCollimatorData_get_record_impacts(el);
        cg->record_exits = TransparentCollimatorData_get_record_exits(el);
    }
    return cg;
}

/*gpufun*/
void TransparentCollimator_free(CollimatorGeometry restrict cg){
    if (cg->side != -1){
        destroy_jaw(cg->segments_L);
    }
    if (cg->side != 1){
        destroy_jaw(cg->segments_R);
    }
    free(cg);
}


/*gpufun*/
void TransparentCollimator_track_local_particle(TransparentCollimatorData el, LocalParticle* part0){
    int8_t active = TransparentCollimatorData_get_active(el);
    active       *= TransparentCollimatorData_get__tracking(el);
    double const length = TransparentCollimatorData_get_length(el);

    // Get geometry
    CollimatorGeometry cg;
    if (active){
        cg = TransparentCollimator_init_geometry(el, part0);
    }

    //start_per_particle_block (part0->part)
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
                int8_t is_hit = hit_jaws_check_and_transform(part, cg);

                // Transform back to the lab frame
                hit_jaws_transform_back(is_hit, part, cg);
                LocalParticle_add_to_s(part, s_coll);
            }
        }
    //end_per_particle_block
    if (active){
        TransparentCollimator_free(cg);
    }
}

#endif /* XCOLL_TRANSPARENT_COLL_H */
