// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_ABSORBER_H
#define XCOLL_ABSORBER_H

/*gpufun*/
CollimatorGeometry BlackAbsorber_init_geometry(BlackAbsorberData el, LocalParticle* part0, int8_t active){
    CollimatorGeometry cg = (CollimatorGeometry) malloc(sizeof(CollimatorGeometry_));
    if (active){ // This is needed in order to avoid that the initialisation is called during a twiss!
        // Collimator jaws (with tilts)
        cg->jaw_LU = BlackAbsorberData_get__jaw_LU(el);
        cg->jaw_LD = BlackAbsorberData_get__jaw_LD(el);
        cg->jaw_RU = BlackAbsorberData_get__jaw_RU(el);
        cg->jaw_RD = BlackAbsorberData_get__jaw_RD(el);
        // TODO: need shortening of active length!
        cg->length = BlackAbsorberData_get_length(el);
        cg->side   = BlackAbsorberData_get__side(el);
        // Get angles of jaws
        cg->sin_zL = BlackAbsorberData_get__sin_zL(el);
        cg->cos_zL = BlackAbsorberData_get__cos_zL(el);
        cg->sin_zR = BlackAbsorberData_get__sin_zR(el);
        cg->cos_zR = BlackAbsorberData_get__cos_zR(el);
        cg->sin_zDiff = BlackAbsorberData_get__sin_zDiff(el);
        cg->cos_zDiff = BlackAbsorberData_get__cos_zDiff(el);
        cg->jaws_parallel = BlackAbsorberData_get__jaws_parallel(el);
        // Impact table
        cg->record = BlackAbsorberData_getp_internal_record(el, part0);
        cg->record_index = NULL;
        cg->record_touches = 0;
        if (cg->record){
            cg->record_index = CollimatorImpactsData_getp__index(cg->record);
            cg->record_touches = BlackAbsorberData_get_record_touches(el);
        }
    }

    return cg;
}


/*gpufun*/
void BlackAbsorber_track_local_particle(BlackAbsorberData el, LocalParticle* part0){

    // Collimator active and length
    int8_t active = BlackAbsorberData_get_active(el);
    active       *= BlackAbsorberData_get__tracking(el);
    double const length = BlackAbsorberData_get_length(el);

    // Get geometry
    CollimatorGeometry cg      = BlackAbsorber_init_geometry(el, part0, active);
    int8_t record_interactions = BlackAbsorberData_get_record_interactions(el);

    //start_per_particle_block (part0->part)
        if (!active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_tracking) {

                // Check if hit on jaws
                int8_t is_hit = hit_jaws_check_and_transform(part, cg);

                if (is_hit != 0){
                    LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);
                    if (record_interactions) {
                        CollimatorImpactsData_log(cg->record, cg->record_index, part, XC_ABSORBED);  // In coll jaw reference frame
                    }
                }

                // Transform back to the lab frame
                hit_jaws_transform_back(is_hit, part, cg);
            }
        }
    //end_per_particle_block

}

#endif /* XCOLL_COLL_GEOM_H */
