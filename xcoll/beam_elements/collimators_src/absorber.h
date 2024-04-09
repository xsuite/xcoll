// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_ABSORBER_H
#define XCOLL_ABSORBER_H


/*gpufun*/
int8_t died_on_jaw(double x, double jaw_L, double jaw_R, int8_t jaws_parallel, LocalParticle* part, double sin_zDiff, double cos_zDiff){
    if (jaws_parallel){
        if (x >= jaw_L){
            return (int8_t) 1;   // Hit on left jaw
        } else if (x <= jaw_R){
            return (int8_t) -1;  // Hit on right jaw
        } else {
            return (int8_t) 0;   // Survived
        }
    } else {
        if (x >= jaw_L){
            return (int8_t) 1;   // Hit on left jaw
        } else {
            // Right jaw is angled w.r.t. to left jaw. Hence, we move to the frame of the right jaw by rotation the difference.
            SRotation_single_particle(part, sin_zDiff, cos_zDiff);
            if (x <= jaw_R){
                return (int8_t) -1;  // Hit on right jaw
            } else {
                // Rotate back.
                SRotation_single_particle(part, -sin_zDiff, cos_zDiff);
                return (int8_t) 0;   // Survived
            }
        }
    }
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
    double sin_zDiff, cos_zDiff;
    int8_t jaws_parallel = 0;
    if (fabs(sin_zL-sin_zR) > 1.e-10 || fabs(cos_zL-cos_zR) > 1.e-10 ){
        jaws_parallel = 1;
        sin_zDiff = sin_zL;
        cos_zDiff = cos_zL;
    } else {
        sin_zDiff = sin_zR*cos_zL - cos_zR*sin_zL;
        cos_zDiff = cos_zL*cos_zR + sin_zL*sin_zR;
    }

    // Impact table
    CollimatorImpactsData record = BlackAbsorberData_getp_internal_record(el, part0);
    RecordIndex record_index = NULL;
    if (record){
        record_index = CollimatorImpactsData_getp__index(record);
    }

    //start_per_particle_block (part0->part)
        if (!is_active){
            // Drift full length
            Drift_single_particle(part, length);

        } else {
            // Check collimator initialisation
            int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);

            int8_t is_dead;
            double x_U, x_D, backtrack_length;

            if (is_tracking) {
                // Rotate to frame of left jaw
                SRotation_single_particle(part, sin_zL, cos_zL);

                // Check if hit on the collimator jaws at the front
                x_U = LocalParticle_get_x(part);
                is_dead = died_on_jaw(x_U, jaw_LU, jaw_RU, jaws_parallel, part, sin_zDiff, cos_zDiff);

                if (is_dead == 1){
                    // Died on left upstream jaw
                    LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);
                    CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW_L);
                    // Rotate back from left jaw
                    SRotation_single_particle(part, -sin_zL, cos_zL);

                } else if (is_dead == -1){
                    // Died on right upstream jaw
                    LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);
                    CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW_R);
                    // Rotate back from right jaw
                    SRotation_single_particle(part, -sin_zR, cos_zR);

                } else {
                    // Continue if the particle didn't hit the collimator
                    Drift_single_particle(part, length);

                    // Check if hit on the collimator jaw at the back
                    x_D = LocalParticle_get_x(part);
                    is_dead = died_on_jaw(x_D, jaw_LD, jaw_RD, jaws_parallel, part, sin_zDiff, cos_zDiff);

                    // If dead, backtrack to the particle position of impact
                    if (is_dead == 1){
                        // Died on left downstream jaw
                        backtrack_length = (jaw_LD - x_D) / (jaw_LD - jaw_LU - x_D + x_U) * length;
                        Drift_single_particle(part, -backtrack_length);
                        LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);
                        CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW_L);
                        // Rotate back from left jaw
                        SRotation_single_particle(part, -sin_zL, cos_zL);

                    } else if (is_dead == -1){
                        // Died on right downstream jaw
                        backtrack_length = (jaw_RD - x_D) / (jaw_RD - jaw_RU - x_D + x_U) * length;
                        Drift_single_particle(part, -backtrack_length);
                        LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);
                        CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW_R);
                        // Rotate back from right jaw
                        SRotation_single_particle(part, -sin_zR, cos_zR);

                    } else {
                        // Survived, rotate back from left jaw
                        SRotation_single_particle(part, -sin_zL, cos_zL);
                    }
                }
            }
        }
    //end_per_particle_block

}

#endif
