#ifndef XCOL_COLLIMATOR_H
#define XCOL_COLLIMATOR_H

/*gpufun*/
int64_t is_within_aperture(
                LocalParticle* part,
                double jaw_L, double jaw_R, double jaw_D, double jaw_U){

    double const x = LocalParticle_get_x(part);
    double const y = LocalParticle_get_y(part);
    return (int64_t)((x > jaw_R) && (x < jaw_L) &&
                     (y > jaw_D) && (y < jaw_U) );
}

/*gpufun*/
void drift_for_collimator(LocalParticle* part, double const length){

    double const rpp    = LocalParticle_get_rpp(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = LocalParticle_get_rvv(part) -
                            ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );

}

/*gpufun*/
void rotation_for_collimator(LocalParticle* part,
                             double const sin_z, double const cos_z){

    	double const x  = LocalParticle_get_x(part);
    	double const y  = LocalParticle_get_y(part);
    	double const px = LocalParticle_get_px(part);
    	double const py = LocalParticle_get_py(part);

    	double const x_hat  =  cos_z * x  + sin_z * y;
    	double const y_hat  = -sin_z * x  + cos_z * y;

    	double const px_hat =  cos_z * px + sin_z * py;
    	double const py_hat = -sin_z * px + cos_z * py;


    	LocalParticle_set_x(part, x_hat);
    	LocalParticle_set_y(part, y_hat);

    	LocalParticle_set_px(part, px_hat);
    	LocalParticle_set_py(part, py_hat);

}

void Collimator_track_local_particle(CollimatorData el, LocalParticle* part0){

    double const inactive_front = CollimatorData_get_inactive_front(el);
    double const inactive_back = CollimatorData_get_inactive_back(el);
    double const active_length = CollimatorData_get_active_length(el);

    double const jaw_L = CollimatorData_get_jaw_L(el);
    double const jaw_R = CollimatorData_get_jaw_R(el);
    double const jaw_U = CollimatorData_get_jaw_U(el);
    double const jaw_D = CollimatorData_get_jaw_D(el);

    double const sin_z = CollimatorData_get_sin_z(el);
    double const cos_z = CollimatorData_get_cos_z(el);
    double const dx = CollimatorData_get_dx(el);
    double const dy = CollimatorData_get_dy(el);

    //start_per_particle_block (part0->part)

        int64_t is_alive;

        // Go to collimator reference system
        LocalParticle_add_to_x(part, -dx);
        LocalParticle_add_to_y(part, -dy);
        rotation_for_collimator(part, sin_z, cos_z);

        // Drift before active length
        drift_for_collimator(part, inactive_front);

        // Drifts and checks in the active part
        is_alive = is_within_aperture(part, jaw_L, jaw_R, jaw_D, jaw_U);
        if (is_alive){
            drift_for_collimator(part, active_length);
        }
        is_alive = is_within_aperture(part, jaw_L, jaw_R, jaw_D, jaw_U);

        // Drift after active length
        if (is_alive){
            drift_for_collimator(part, inactive_back);
        }

        // Back from collimator reference system
        rotation_for_collimator(part, -sin_z, cos_z);
        LocalParticle_add_to_x(part, dx);
        LocalParticle_add_to_y(part, dy);

        if (!is_alive){
           LocalParticle_set_state(part, -333);
	    }



    //end_per_particle_block

}

#endif