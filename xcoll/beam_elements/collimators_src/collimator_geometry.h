// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_H
#define XCOLL_COLL_GEOM_H
#include <math.h>
#include <stdio.h>

typedef struct CollimatorGeometry_ {
    // Collimator jaws (with tilts)
    double jaw_LU;
    double jaw_LD;
    double jaw_RU;
    double jaw_RD;
    // TODO: need shortening of active length!
    double length;
    int8_t side  ;
    // Get angles of jaws
    double sin_zL;
    double cos_zL;
    double sin_zR;
    double cos_zR;
    double sin_zDiff;
    double cos_zDiff;
    int8_t jaws_parallel;
    // Impact table
    InteractionRecordData record;
    RecordIndex record_index;
    int8_t record_touches;
} CollimatorGeometry_;
typedef CollimatorGeometry_ *CollimatorGeometry;


// This function expects the particles to be rotated (in the XY plane) to the reference frame of the left jaw
/*gpufun*/
int8_t _hit_on_one_jaw(double x, double jaw_L, double jaw_R, LocalParticle* part, CollimatorGeometry cg){
    if (cg->jaws_parallel){
        if (x >= jaw_L && cg->side != -1){
            return (int8_t) 1;   // Hit on left jaw
        } else if (x <= jaw_R && cg->side != 1){
            return (int8_t) -1;  // Hit on right jaw
        } else {
            return (int8_t) 0;   // Survived
        }

    } else {
        if (x >= jaw_L && cg->side != -1){
            return (int8_t) 1;   // Hit on left jaw
        } else if (cg->side == 1){   // Left-sided collimator
            return (int8_t) 0;   // Survived
        } else {
            // Right jaw is angled w.r.t. to left jaw. Hence, we move to the frame of the right jaw by rotating the difference.
            SRotation_single_particle(part, cg->sin_zDiff, cg->cos_zDiff);
            if (x <= jaw_R){
                return (int8_t) -1;  // Hit on right jaw
            } else {
                // Rotate back to left frame.
                SRotation_single_particle(part, -cg->sin_zDiff, cg->cos_zDiff);
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
int8_t hit_jaws_check_and_transform(LocalParticle* part, CollimatorGeometry cg){
    double x_U, x_D, backtrack_length;

    // Rotate to the frame of the left jaw
    SRotation_single_particle(part, cg->sin_zL, cg->cos_zL);

    // Check if hit on the collimator jaws at the front
    x_U = LocalParticle_get_x(part);
    int8_t is_hit = _hit_on_one_jaw(x_U, cg->jaw_LU, cg->jaw_RU, part, cg);

    if (is_hit == 1){
        // Hit on left upstream jaw
        if (cg->record_touches){
            InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_L);
        }

    } else if (is_hit == -1){
        // Hit on right upstream jaw
        if (cg->record_touches){
            InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_R);
        }

    }

    // Continue if the particle didn't hit the collimator front
    if (is_hit == 0){
        Drift_single_particle(part, cg->length);

        // Check if hit on the collimator jaw at the back
        x_D = LocalParticle_get_x(part);
        is_hit = _hit_on_one_jaw(x_D, cg->jaw_LD, cg->jaw_RD, part, cg);

        // If hit, backtrack to the particle position of impact
        if (is_hit == 1){
            // Hit on left downstream jaw
            backtrack_length = (cg->jaw_LD - x_D) / (cg->jaw_LD - cg->jaw_LU - x_D + x_U) * cg->length;
            Drift_single_particle(part, -backtrack_length);
            if (cg->record_touches){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_L);
            }

        } else if (is_hit == -1){
            // Hit on right downstream jaw
            backtrack_length = (cg->jaw_RD - x_D) / (cg->jaw_RD - cg->jaw_RU - x_D + x_U) * cg->length;
            Drift_single_particle(part, -backtrack_length);
            if (cg->record_touches){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_R);
            }

        } else {
            // If still not hit, return to lab frame
            SRotation_single_particle(part, -cg->sin_zL, cg->cos_zL);

        }
    }

    return is_hit;
}


/*gpufun*/
void hit_jaws_transform_back(int8_t is_hit, LocalParticle* part, CollimatorGeometry cg){
    if (is_hit == 1){
        // Hit on left jaw
        SRotation_single_particle(part, -cg->sin_zL, cg->cos_zL);
    } else if (is_hit == -1){
        // Hit on right jaw
        SRotation_single_particle(part, -cg->sin_zR, cg->cos_zR);
    }
}


#endif /* XCOLL_COLL_GEOM_H */
