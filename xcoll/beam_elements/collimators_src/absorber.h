// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_ABSORBER_H
#define XCOLL_ABSORBER_H


/*gpufun*/
void BlackAbsorber_track_local_particle(BlackAbsorberData el, LocalParticle* part0){

    // Collimator active and length
    int8_t is_active = BlackAbsorberData_get_active(el);
    is_active       *= BlackAbsorberData_get__tracking(el);
    double const length = BlackAbsorberData_get_length(el);

    // Impact table
    RecordIndex record_index = NULL;
    int8_t record_touches = 0;
    int8_t record_interactions = 0;
    CollimatorImpactsData record = BlackAbsorberData_getp_internal_record(el, part0);
    if (record){
        record_index = CollimatorImpactsData_getp__index(record);
        record_touches = BlackAbsorberData_get_record_touches(el);
        record_interactions = BlackAbsorberData_get_record_interactions(el);
    }

    //start_per_particle_block (part0->part)
        if (!is_active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);

            if (is_tracking) {

                // Check if hit on jaws
                int8_t is_hit = hit_jaws_check_and_transform(part, (BaseCollimatorData) el, record, record_index, record_touches);

                if (is_hit == 1){
                    // Died on left jaw
                    LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);
                    if (record_interactions) {
                        CollimatorImpactsData_log(record, record_index, part, XC_ABSORBED);  // In coll jaw reference frame
                    }

                } else if (is_hit == -1){
                    // Died on right jaw
                    LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);
                    if (record_interactions){
                        CollimatorImpactsData_log(record, record_index, part, XC_ABSORBED);  // In coll jaw reference frame
                    }
                }

                // Transform back to the lab frame
                hit_jaws_transform_back(is_hit, part, (BaseCollimatorData) el);
            }
        }
    //end_per_particle_block

}

#endif /* XCOLL_COLL_GEOM_H */
