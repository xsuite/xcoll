#ifndef XCOLL_ABSORBER_H
#define XCOLL_ABSORBER_H

/*gpufun*/
int64_t is_within_aperture(LocalParticle* part,
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

void BlackAbsorber_track_local_particle(BlackAbsorberData el, LocalParticle* part0){

    // Collimator active and length
    int8_t const is_active = BlackAbsorberData_get_active(el);
    double const inactive_front = BlackAbsorberData_get_inactive_front(el);
    double const inactive_back = BlackAbsorberData_get_inactive_back(el);
    double const active_length = BlackAbsorberData_get_active_length(el);
    // Collimator jaws
    double const jaw_F_L = BlackAbsorberData_get_jaw_F_L(el);
    double const jaw_F_R = BlackAbsorberData_get_jaw_F_R(el);
    double const jaw_B_L = BlackAbsorberData_get_jaw_B_L(el);
    double const jaw_B_R = BlackAbsorberData_get_jaw_B_R(el);
    double const jaw_U = BlackAbsorberData_get_jaw_U(el);
    double const jaw_D = BlackAbsorberData_get_jaw_D(el);
    // Collimator reference frame
    double const sin_z = BlackAbsorberData_get_sin_z(el);
    double const cos_z = BlackAbsorberData_get_cos_z(el);
    double const dx = BlackAbsorberData_get_dx(el);
    double const dy = BlackAbsorberData_get_dy(el);
    // Impact table
    int8_t  const record_impacts = BlackAbsorberData_get__record_impacts(el);
    if (record_impacts){
        int64_t const capacity = BlackAbsorberData_get_impacts__capacity(el);
        assert(record_index < capacity - 1);
    }


    //start_per_particle_block (part0->part)

        // Go to collimator reference system (centered around orbit)
        LocalParticle_add_to_x(part, -dx);
        LocalParticle_add_to_y(part, -dy);
        rotation_for_collimator(part, sin_z, cos_z);

        if (!is_active){

            // If collimator not active, replace with drift
            drift_for_collimator(part, inactive_front + active_length + inactive_back);

        } else {

            int64_t is_alive;

            // Store transversal coordinates for potential backtracking later
            double x_F = LocalParticle_get_x(part);
            double y_F = LocalParticle_get_y(part);
            
            // Drift inactive length before jaw
            drift_for_collimator(part, inactive_front);

            // Check if hit on the collimator jaw at the front
            is_alive = is_within_aperture(part, jaw_F_L, jaw_F_R, jaw_D, jaw_U);

            // Continue if the particle didn't hit the collimator
            if (is_alive){

                // Drift the jaw length
                drift_for_collimator(part, active_length);

                // Check if hit on the collimator jaw at the back
                is_alive = is_within_aperture(part, jaw_B_L, jaw_B_R, jaw_D, jaw_U);

                // TODO: is there a performance difference with nesting the ifs or not?
                // Continue if the particle didn't hit the collimator
                if (is_alive){
                    
                    // Drift inactive length after jaw
                    drift_for_collimator(part, inactive_back);
                    
                } else {
                    
                    // Backtrack to the particle position of impact
                    // This should only be done if the particle did NOT hit the front jaw
                    double x_B = LocalParticle_get_x(part);
                    double y_B = LocalParticle_get_y(part);
                    double length;
                    
                    if (x_B > 0){        // Left jaw
                        length = (jaw_B_L - x_B) / (jaw_B_L - jaw_F_L - x_B + x_F) * active_length;
                    } else if (x_B < 0){ // Right jaw
                        length = (jaw_B_R - x_B) / (jaw_B_R - jaw_F_R - x_B + x_F) * active_length;
                    } else if (y_B > 0){ // Upper jaw
                        length = (y_B - jaw_U) / (y_B - y_F) * active_length;
                    } else if (y_B < 0){ // Lower jaw
                        length = (y_B - jaw_D) / (y_B - y_F) * active_length;
                    } else {
                        length = 0;
                    }

                    drift_for_collimator(part, length);

                }
            }

            if (!is_alive){

                LocalParticle_set_state(part, -333);
                // Record data
                if (record_impacts){
//                     // TODO: GPU-proof way (though did not work):
//                     BlackAbsorberData_add_to_impacts__row_id(el, 1);
//                     int64_t record_index = BlackAbsorberData_get_impacts__row_id(el);

                    int64_t record_index = BlackAbsorberData_get_impacts__row_id(el);
                    record_index = record_index + 1;
                    BlackAbsorberData_set_impacts__row_id(el, record_index);

                    double mass_ratio = LocalParticle_get_charge_ratio(part) / LocalParticle_get_chi(part);
                    double m = LocalParticle_get_mass0(part) * mass_ratio;
                    double q = LocalParticle_get_charge_ratio(part) * LocalParticle_get_q0(part);
                    double energy = (LocalParticle_get_ptau(part) + 1) * mass_ratio * LocalParticle_get_p0c(part);

                    BlackAbsorberData_set_impacts_at_element(el, record_index, LocalParticle_get_at_element(part));
                    BlackAbsorberData_set_impacts_s(el, record_index, LocalParticle_get_s(part));
                    // TODO: turn is zero; does that make sense?
                    BlackAbsorberData_set_impacts_turn(el, record_index, LocalParticle_get_at_turn(part));
//                     BlackAbsorberData_set_impacts_interaction_type(el, record_index, 0);

                    BlackAbsorberData_set_impacts_id_parent(el, record_index, LocalParticle_get_particle_id(part));
                    BlackAbsorberData_set_impacts_id_child(el, record_index, -1);
                    BlackAbsorberData_set_impacts_x_parent(el, record_index, LocalParticle_get_x(part));
                    BlackAbsorberData_set_impacts_px_parent(el, record_index, LocalParticle_get_px(part));
                    BlackAbsorberData_set_impacts_y_parent(el, record_index, LocalParticle_get_y(part));
                    BlackAbsorberData_set_impacts_py_parent(el, record_index, LocalParticle_get_py(part));
                    BlackAbsorberData_set_impacts_zeta_parent(el, record_index, LocalParticle_get_zeta(part));
                    BlackAbsorberData_set_impacts_delta_parent(el, record_index, LocalParticle_get_delta(part));
                    BlackAbsorberData_set_impacts_energy_parent(el, record_index, energy);
                    BlackAbsorberData_set_impacts_m_parent(el, record_index, m);
                    BlackAbsorberData_set_impacts_q_parent(el, record_index, q);
//                     BlackAbsorberData_set_impacts_a_parent(el, record_index, LocalParticle_get_a(part));
//                     BlackAbsorberData_set_impacts_pdgid_parent(el, record_index, LocalParticle_get_pdgid(part));
                }
            }

        }

        // Move back from collimator reference system
        rotation_for_collimator(part, -sin_z, cos_z);
        LocalParticle_add_to_x(part, dx);
        LocalParticle_add_to_y(part, dy);

    //end_per_particle_block

}

#endif