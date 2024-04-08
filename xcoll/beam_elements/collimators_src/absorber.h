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
    double const length = BlackAbsorberData_get_length(el);

    // Collimator jaws
    double const jaw_L  = BlackAbsorberData_get_jaw_L(el);
    double const jaw_R  = BlackAbsorberData_get_jaw_R(el);
    double const sin_yL = BlackAbsorberData_get_sin_yL(el);
    double const sin_yR = BlackAbsorberData_get_sin_yR(el);
    double const jaw_LU = jaw_L - sin_yL*length/2.;
    double const jaw_LD = jaw_L + sin_yL*length/2.;
    double const jaw_RU = jaw_R - sin_yR*length/2.;
    double const jaw_RD = jaw_R + sin_yR*length/2.;
    // TODO: need shortening of active length!
    // Collimator reference frame
    double const sin_zL = BlackAbsorberData_get_sin_zL(el);
    double const cos_zL = BlackAbsorberData_get_cos_zL(el);
    double const sin_zR = BlackAbsorberData_get_sin_zR(el);
    double const cos_zR = BlackAbsorberData_get_cos_zR(el);
    if (fabs(sin_zL-sin_zR) > 1.e-10 || fabs(cos_zL-cos_zR) > 1.e-10 ){
        kill_all_particles(part0, XC_ERR_NOT_IMPLEMENTED);
    };

    // Impact table
    CollimatorImpactsData record = BlackAbsorberData_getp_internal_record(el, part0);
    RecordIndex record_index = NULL;
    if (record){
        record_index = CollimatorImpactsData_getp__index(record);
    }

    //start_per_particle_block (part0->part)
    
        double s_collimator = LocalParticle_get_s(part);

        // Go to collimator reference system (centered around orbit)
        SRotation_single_particle(part, sin_zL, cos_zL);

        int64_t is_alive = 1;

        if (!is_active){

            // If collimator not active, replace with drift
            Drift_single_particle(part, length);

        } else {

            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
            if (is_tracking) {

            // Store transversal coordinates for potential backtracking later
            double x_F = LocalParticle_get_x(part);
//             double y_F = LocalParticle_get_y(part);

            // Check if hit on the collimator jaw at the front
            is_alive = is_within_aperture(part, jaw_LU, jaw_RU);

            // Continue if the particle didn't hit the collimator
            if (is_alive){

                // Drift the jaw length
                Drift_single_particle(part, length);

                // Check if hit on the collimator jaw at the back
                is_alive = is_within_aperture(part, jaw_LD, jaw_RD);

                // TODO: is there a performance difference with nesting the ifs or not?
                // Continue if the particle didn't hit the collimator
                if (!is_alive){

                    // Backtrack to the particle position of impact
                    // This should only be done if the particle did NOT hit the front jaw
                    double x_B = LocalParticle_get_x(part);
//                     double y_B = LocalParticle_get_y(part);
                    double backtrack_length;

                    if (x_B > 0){        // Left jaw
                        backtrack_length = (jaw_LD - x_B) / (jaw_LD - jaw_LU - x_B + x_F) * length;
                    } else if (x_B < 0){ // Right jaw
                        backtrack_length = (jaw_RD - x_B) / (jaw_RD - jaw_RU - x_B + x_F) * length;
                    // TODO: check this
//                     } else if (y_B > 0){ // Upper jaw
//                         length = (y_B - jaw_U) / (y_B - y_F) * length;
//                     } else if (y_B < 0){ // Lower jaw
//                         length = (y_B - jaw_D) / (y_B - y_F) * length;
                    } else {
                        backtrack_length = 0;
                    }
                    Drift_single_particle(part, -backtrack_length);

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

    //end_per_particle_block

}

#endif
