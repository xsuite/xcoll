// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_H
#define XCOLL_COLL_GEOM_H
#include <math.h>
#include <stdio.h>


// This function expects the particles to be rotated (in the XY plane) to the reference frame of the left jaw
/*gpufun*/
int8_t _hit_on_one_jaw(double x, int8_t side, double jaw_L, double jaw_R, int8_t jaws_parallel,
                      LocalParticle* part, BaseCollimatorData el){
    if (jaws_parallel){
        if (x >= jaw_L && side != -1){
            return (int8_t) 1;   // Hit on left jaw
        } else if (x <= jaw_R && side != 1){
            return (int8_t) -1;  // Hit on right jaw
        } else {
            return (int8_t) 0;   // Survived
        }

    } else {
        if (x >= jaw_L && side != -1){
            return (int8_t) 1;   // Hit on left jaw
        } else if (side == 1){   // Left-sided collimator
            return (int8_t) 0;   // Survived
        } else {
            // Right jaw is angled w.r.t. to left jaw. Hence, we move to the frame of the right jaw by rotating the difference.
            double const sin_zDiff = BaseCollimatorData_get__sin_zDiff(el);
            double const cos_zDiff = BaseCollimatorData_get__cos_zDiff(el);
            SRotation_single_particle(part, sin_zDiff, cos_zDiff);
            if (x <= jaw_R){
                return (int8_t) -1;  // Hit on right jaw
            } else {
                // Rotate back to left frame.
                SRotation_single_particle(part, -sin_zDiff, cos_zDiff);
                return (int8_t) 0;   // Survived
            }
        }
    }
}



// This function checks if a particle hits a jaw (and which).
// Return value: 0 (no hit), 1 (hit on left jaw), -1 (hit on right jaw).
// Furthermore, the particle is moved to the location where it hits the jaw (drifted to the end if no hit),
//              and transformed to the reference frame of that jaw.

/*gpufun*/
int8_t hit_jaws_check_and_transform(LocalParticle* part, BaseCollimatorData el, CollimatorImpactsData record,
                                    RecordIndex record_index, int8_t record_touches){
    // Get data
    // ========

    // Collimator jaws (with tilts)
    double const jaw_LU = BaseCollimatorData_get__jaw_LU(el);
    double const jaw_LD = BaseCollimatorData_get__jaw_LD(el);
    double const jaw_RU = BaseCollimatorData_get__jaw_RU(el);
    double const jaw_RD = BaseCollimatorData_get__jaw_RD(el);
    // TODO: need shortening of active length!
    double const length = BaseCollimatorData_get_length(el);
    int8_t const side   = BaseCollimatorData_get__side(el);

    // Get angles of jaws
    double const sin_zL = BaseCollimatorData_get__sin_zL(el);
    double const cos_zL = BaseCollimatorData_get__cos_zL(el);
    int8_t const jaws_parallel = BaseCollimatorData_get__jaws_parallel(el);

    // Do the tracking
    // ==============

    double x_U, x_D, backtrack_length;

    // Rotate to the frame of the left jaw
    SRotation_single_particle(part, sin_zL, cos_zL);

    // Check if hit on the collimator jaws at the front
    x_U = LocalParticle_get_x(part);
    int8_t is_hit = _hit_on_one_jaw(x_U, side, jaw_LU, jaw_RU, jaws_parallel, part, el);

    // Log the hit
    if (record_touches){
        if (is_hit == 1){
            // Hit on left upstream jaw
            CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW_L);
    
        } else if (is_hit == -1){
            // Hit on right upstream jaw
            CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW_R);
    
        }
    }

    // Continue if the particle didn't hit the collimator front
    if (is_hit == 0){
        Drift_single_particle(part, length);

        // Check if hit on the collimator jaw at the back
        x_D = LocalParticle_get_x(part);
        is_hit = _hit_on_one_jaw(x_D, side, jaw_LD, jaw_RD, jaws_parallel, part, el);

        // If hit, backtrack to the particle position of impact
        if (is_hit == 1){
            // Hit on left downstream jaw
            backtrack_length = (jaw_LD - x_D) / (jaw_LD - jaw_LU - x_D + x_U) * length;
            Drift_single_particle(part, -backtrack_length);
            if (record_touches){
                CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW_L);
            }

        } else if (is_hit == -1){
            // Hit on right downstream jaw
            backtrack_length = (jaw_RD - x_D) / (jaw_RD - jaw_RU - x_D + x_U) * length;
            Drift_single_particle(part, -backtrack_length);
            if (record_touches){
                CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW_R);
            }

        } else {
            // If still not hit, return to lab frame
            SRotation_single_particle(part, -sin_zL, cos_zL);

        }
    }

    return is_hit;
}


/*gpufun*/
void hit_jaws_transform_back(int8_t is_hit, LocalParticle* part, BaseCollimatorData el){
    if (is_hit == 1){
        // Hit on left jaw
        double const sin_zL = BaseCollimatorData_get__sin_zL(el);
        double const cos_zL = BaseCollimatorData_get__cos_zL(el);
        SRotation_single_particle(part, -sin_zL, cos_zL);
    } else if (is_hit == -1){
        // Hit on right jaw
        double const sin_zR = BaseCollimatorData_get__sin_zR(el);
        double const cos_zR = BaseCollimatorData_get__cos_zR(el);
        SRotation_single_particle(part, -sin_zR, cos_zR);
    }
}


#endif /* XCOLL_COLL_GEOM_H */
