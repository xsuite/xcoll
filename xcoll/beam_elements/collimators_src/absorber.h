// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_ABSORBER_H
#define XCOLL_ABSORBER_H

/*gpufun*/
int64_t is_within_aperture(LocalParticle* part, double jaw_L, double jaw_R){
    double const x = LocalParticle_get_x(part);
    return (int64_t)((x > jaw_R) && (x < jaw_L) );
}

/*gpufun*/
void BlackAbsorber_track_local_particle(BlackAbsorberData el, LocalParticle* part0){

    // Collimator active and length
    int8_t is_active = BlackAbsorberData_get_active(el);
    is_active       *= BlackAbsorberData_get__tracking(el);
    double const inactive_front = BlackAbsorberData_get_inactive_front(el);
    double const inactive_back = BlackAbsorberData_get_inactive_back(el);
    double const active_length = BlackAbsorberData_get_active_length(el);
    // Collimator jaws
    double const jaw_L  = BlackAbsorberData_get_jaw_L(el);
    double const jaw_R  = BlackAbsorberData_get_jaw_R(el);
    double const sin_yL = BlackAbsorberData_get_sin_yL(el);
    double const sin_yR = BlackAbsorberData_get_sin_yR(el);
    double const jaw_LU = jaw_L - sin_yL*active_length/2.;
    double const jaw_LD = jaw_L + sin_yL*active_length/2.;
    double const jaw_RU = jaw_R - sin_yR*active_length/2.;
    double const jaw_RD = jaw_R + sin_yR*active_length/2.;
    // TODO: need shortening of active length!
    // Collimator reference frame
    double const sin_zL = BlackAbsorberData_get_sin_zL(el);
    double const cos_zL = BlackAbsorberData_get_cos_zL(el);
    double const sin_zR = BlackAbsorberData_get_sin_zR(el);
    double const cos_zR = BlackAbsorberData_get_cos_zR(el);
    if (fabs(sin_zL-sin_zR) > 1.e-10 || fabs(cos_zL-cos_zR) > 1.e-10 ){
        kill_all_particles(part0, XC_ERR_NOT_IMPLEMENTED);
    };
    double const dx = BlackAbsorberData_get_ref_x(el);
    double const dy = BlackAbsorberData_get_ref_y(el);
    // Impact table
    CollimatorImpactsData record = BlackAbsorberData_getp_internal_record(el, part0);
    RecordIndex record_index = NULL;
    if (record){
        record_index = CollimatorImpactsData_getp__index(record);
    }

    //start_per_particle_block (part0->part)
    
        double s_collimator = LocalParticle_get_s(part);

        // Go to collimator reference system (centered around orbit)
        XYShift_single_particle(part, dx, dy);
        SRotation_single_particle(part, sin_zL, cos_zL);

        int64_t is_alive = 1;

        if (!is_active){

            // If collimator not active, replace with drift
            Drift_single_particle(part, inactive_front+active_length+inactive_back);

        } else {

            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
            if (is_tracking) {
           
            // Drift inactive length before jaw
            Drift_single_particle(part, inactive_front);

            // Store transversal coordinates for potential backtracking later
            double x_F = LocalParticle_get_x(part);
//             double y_F = LocalParticle_get_y(part);

            // Check if hit on the collimator jaw at the front
            is_alive = is_within_aperture(part, jaw_LU, jaw_RU);

            // Continue if the particle didn't hit the collimator
            if (is_alive){

                // Drift the jaw length
                Drift_single_particle(part, active_length);

                // Check if hit on the collimator jaw at the back
                is_alive = is_within_aperture(part, jaw_LD, jaw_RD);

                // TODO: is there a performance difference with nesting the ifs or not?
                // Continue if the particle didn't hit the collimator
                if (is_alive){

                    // Drift inactive length after jaw
                    Drift_single_particle(part, inactive_back);

                } else {

                    // Backtrack to the particle position of impact
                    // This should only be done if the particle did NOT hit the front jaw
                    double x_B = LocalParticle_get_x(part);
//                     double y_B = LocalParticle_get_y(part);
                    double length;

                    if (x_B > 0){        // Left jaw
                        length = (jaw_LD - x_B) / (jaw_LD - jaw_LU - x_B + x_F) * active_length;
                    } else if (x_B < 0){ // Right jaw
                        length = (jaw_RD - x_B) / (jaw_RD - jaw_RU - x_B + x_F) * active_length;
                    // TODO: check this
//                     } else if (y_B > 0){ // Upper jaw
//                         length = (y_B - jaw_U) / (y_B - y_F) * active_length;
//                     } else if (y_B < 0){ // Lower jaw
//                         length = (y_B - jaw_D) / (y_B - y_F) * active_length;
                    } else {
                        length = 0;
                    }
                    Drift_single_particle(part, -length);

                }
            }
            }
        }

        // Update dead particles
        if (!is_alive){

            LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);

            // Record impact data
            CollimatorImpactsData_log(record, record_index, part, XC_ABSORBED);
        }

        // Move back from collimator reference system
        SRotation_single_particle(part, -sin_zL, cos_zL);
        XYShift_single_particle(part, -dx, -dy);

    //end_per_particle_block

}

#endif
