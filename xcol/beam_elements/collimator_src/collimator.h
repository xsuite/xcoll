#ifndef XCOL_COLLIMATOR_H
#define XCOL_COLLIMATOR_H

/*gpufun*/
int64_t is_within_aperture(
                LocalParticle* part,
                double a_min, double a_max, double b_min, double b_max){

    double const x = LocalParticle_get_x(part);
    double const y = LocalParticle_get_y(part);
    return (int64_t)((x >= a_min) && (x <= a_max) &&
                     (y >= b_min) && (y <= b_max) );
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

    double const inactive_length_at_start =
                        CollimatorData_get_inactive_length_at_start(el);
    double const inactive_length_at_end =
                        CollimatorData_get_inactive_length_at_end(el);
    int64_t const n_slices = CollimatorData_get_n_slices(el);

    double const slice_length = CollimatorData_get_active_length(el)/n_slices;

    double const a_max = CollimatorData_get_a_max(el);
    double const a_min = CollimatorData_get_a_min(el);
    double const b_max = CollimatorData_get_b_max(el);
    double const b_min = CollimatorData_get_b_min(el);

    double const sin_z = CollimatorData_get_sin_z(el);
    double const cos_z = CollimatorData_get_cos_z(el);

    //start_per_particle_block (part0->part)

        int64_t is_alive;

        // Drift before active length

        rotation_for_collimator(part, sin_z, cos_z);

        drift_for_collimator(part, inactive_length_at_start);

        for (int64_t islice=0; islice<n_slices; islice++){

            is_alive = is_within_aperture(part, a_min, a_max, b_min, b_max);
            if (!is_alive){break;}

            drift_for_collimator(part, slice_length);

        }
        is_alive = is_within_aperture(part, a_min, a_max, b_min, b_max);

        if (is_alive){
            drift_for_collimator(part, inactive_length_at_end);
        }

        rotation_for_collimator(part, -sin_z, cos_z);

        if (!is_alive){
           LocalParticle_set_state(part, -333);
	    }



    //end_per_particle_block

}

#endif