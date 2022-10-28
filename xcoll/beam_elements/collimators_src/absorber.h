#ifndef XCOLL_ABSORBER_H
#define XCOLL_ABSORBER_H

/*gpufun*/
int64_t is_within_aperture(LocalParticle* part, double jaw_L, double jaw_R){

    double const x = LocalParticle_get_x(part);
    return (int64_t)((x > jaw_R) && (x < jaw_L) );
}

/*gpufun*/
void drift_for_collimator(LocalParticle* part, double const length){

    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

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

/*gpufun*/
void BlackAbsorber_track_local_particle(BlackAbsorberData el, LocalParticle* part0){

    // Collimator active and length
    int8_t const is_active = BlackAbsorberData_get__active(el);
    double const inactive_front = BlackAbsorberData_get_inactive_front(el);
    double const inactive_back = BlackAbsorberData_get_inactive_back(el);
    double const active_length = BlackAbsorberData_get_active_length(el);
    // Collimator jaws
    double const jaw_F_L = BlackAbsorberData_get_jaw_F_L(el);
    double const jaw_F_R = BlackAbsorberData_get_jaw_F_R(el);
    double const jaw_B_L = BlackAbsorberData_get_jaw_B_L(el);
    double const jaw_B_R = BlackAbsorberData_get_jaw_B_R(el);
    // Collimator reference frame
    double const sin_z = BlackAbsorberData_get_sin_z(el);
    double const cos_z = BlackAbsorberData_get_cos_z(el);
    double const dx = BlackAbsorberData_get_dx(el);
    double const dy = BlackAbsorberData_get_dy(el);
    // Impact table
    CollimatorImpactsData record = BlackAbsorberData_getp_internal_record(el, part0);
    RecordIndex record_index = NULL;
    if (record){
        record_index = CollimatorImpactsData_getp__index(record);
    }

    //start_per_particle_block (part0->part)

        // Go to collimator reference system (centered around orbit)
        LocalParticle_add_to_x(part, -dx);
        LocalParticle_add_to_y(part, -dy);
        rotation_for_collimator(part, sin_z, cos_z);

        int64_t is_alive = 1;

        if (!is_active){

            // If collimator not active, replace with drift
            drift_for_collimator(part, inactive_front + active_length + inactive_back);

        } else {
           
            // Drift inactive length before jaw
            drift_for_collimator(part, inactive_front);

            // Store transversal coordinates for potential backtracking later
            double x_F = LocalParticle_get_x(part);
//             double y_F = LocalParticle_get_y(part);

            // Check if hit on the collimator jaw at the front
            is_alive = is_within_aperture(part, jaw_F_L, jaw_F_R);

            // Continue if the particle didn't hit the collimator
            if (is_alive){

                // Drift the jaw length
                drift_for_collimator(part, active_length);

                // Check if hit on the collimator jaw at the back
                is_alive = is_within_aperture(part, jaw_B_L, jaw_B_R);

                // TODO: is there a performance difference with nesting the ifs or not?
                // Continue if the particle didn't hit the collimator
                if (is_alive){
                    
                    // Drift inactive length after jaw
                    drift_for_collimator(part, inactive_back);
                    
                } else {
                    
                    // Backtrack to the particle position of impact
                    // This should only be done if the particle did NOT hit the front jaw
                    double x_B = LocalParticle_get_x(part);
//                     double y_B = LocalParticle_get_y(part);
                    double length;
                    
                    if (x_B > 0){        // Left jaw
                        length = (jaw_B_L - x_B) / (jaw_B_L - jaw_F_L - x_B + x_F) * active_length;
                    } else if (x_B < 0){ // Right jaw
                        length = (jaw_B_R - x_B) / (jaw_B_R - jaw_F_R - x_B + x_F) * active_length;
                    // TODO: check this
//                     } else if (y_B > 0){ // Upper jaw
//                         length = (y_B - jaw_U) / (y_B - y_F) * active_length;
//                     } else if (y_B < 0){ // Lower jaw
//                         length = (y_B - jaw_D) / (y_B - y_F) * active_length;
                    } else {
                        length = 0;
                    }
                    drift_for_collimator(part, -length);

                }
            }
        }

        // Move back from collimator reference system
        rotation_for_collimator(part, -sin_z, cos_z);
        LocalParticle_add_to_x(part, dx);
        LocalParticle_add_to_y(part, dy);

        // Update dead particles
        if (!is_alive){

            LocalParticle_set_state(part, -333);

            // Record impact data
            if (record){
                // Get a slot in the record (this is thread safe)
                int64_t i_slot = RecordIndex_get_slot(record_index);
                // The returned slot id is negative if record is NULL or if record is full

                if (i_slot>=0){
                    double mass_ratio = LocalParticle_get_charge_ratio(part) / LocalParticle_get_chi(part);
                    double energy = ( LocalParticle_get_ptau(part) + 1 / LocalParticle_get_beta0(part)
                                     ) * mass_ratio * LocalParticle_get_p0c(part);

                    CollimatorImpactsData_set_at_element(record, i_slot, LocalParticle_get_at_element(part));
                    CollimatorImpactsData_set_at_turn(record, i_slot, LocalParticle_get_at_turn(part));
                    CollimatorImpactsData_set_s(record, i_slot, LocalParticle_get_s(part));
                    CollimatorImpactsData_set_interaction_type(record, i_slot, -1);

                    CollimatorImpactsData_set_parent_id(record, i_slot, LocalParticle_get_particle_id(part));
                    CollimatorImpactsData_set_parent_x(record, i_slot, LocalParticle_get_x(part));
                    CollimatorImpactsData_set_parent_px(record, i_slot, LocalParticle_get_px(part));
                    CollimatorImpactsData_set_parent_y(record, i_slot, LocalParticle_get_y(part));
                    CollimatorImpactsData_set_parent_py(record, i_slot, LocalParticle_get_py(part));
                    CollimatorImpactsData_set_parent_zeta(record, i_slot, LocalParticle_get_zeta(part));
                    CollimatorImpactsData_set_parent_delta(record, i_slot, LocalParticle_get_delta(part));
                    CollimatorImpactsData_set_parent_energy(record, i_slot, energy);
                    // TODO: particle info
                    CollimatorImpactsData_set_parent_mass(record, i_slot, -1);
                    CollimatorImpactsData_set_parent_charge(record, i_slot, 1);
                    CollimatorImpactsData_set_parent_z(record, i_slot, -1);
                    CollimatorImpactsData_set_parent_a(record, i_slot, -1);
                    CollimatorImpactsData_set_parent_pdgid(record, i_slot, -1);

                    CollimatorImpactsData_set_child_id(record, i_slot, -1);
                    // We need to fill in child data, or the arrays will not have the same length
                    // (when secondaries are be made elsewhere in other collimators)
                    CollimatorImpactsData_set_child_x(record, i_slot, -1);
                    CollimatorImpactsData_set_child_px(record, i_slot, -1);
                    CollimatorImpactsData_set_child_y(record, i_slot, -1);
                    CollimatorImpactsData_set_child_py(record, i_slot, -1);
                    CollimatorImpactsData_set_child_zeta(record, i_slot, -1);
                    CollimatorImpactsData_set_child_delta(record, i_slot, -1);
                    CollimatorImpactsData_set_child_energy(record, i_slot, -1);
                    CollimatorImpactsData_set_child_mass(record, i_slot, -1);
                    CollimatorImpactsData_set_child_charge(record, i_slot, -1);
                    CollimatorImpactsData_set_child_z(record, i_slot, -1);
                    CollimatorImpactsData_set_child_a(record, i_slot, -1);
                    CollimatorImpactsData_set_child_pdgid(record, i_slot, -1);
                }
            }
        }

    //end_per_particle_block

}

#endif
