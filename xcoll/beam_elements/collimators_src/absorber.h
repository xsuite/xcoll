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
    double const jaw_LU = BlackAbsorberData_get_jaw_LU(el);
    double const jaw_RU = BlackAbsorberData_get_jaw_RU(el);
    double const jaw_LD = BlackAbsorberData_get_jaw_LD(el);
    double const jaw_RD = BlackAbsorberData_get_jaw_RD(el);
    // Collimator reference frame
    double const sin_z = BlackAbsorberData_get_sin_zL(el);
    double const cos_z = BlackAbsorberData_get_cos_zL(el);
    double const dx = BlackAbsorberData_get_ref_x(el);
    double const dy = BlackAbsorberData_get_ref_y(el);
    // Impact table
    CollimatorImpactsData record = BlackAbsorberData_getp_internal_record(el, part0);
    RecordIndex record_index = NULL;
    if (record){
        record_index = CollimatorImpactsData_getp__index(record);
    }

    //start_per_particle_block (part0->part)

        // Go to collimator reference system (centered around orbit)
        XYShift_single_particle(part, dx, dy);
        SRotation_single_particle(part, sin_z, cos_z);

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

        // Move back from collimator reference system
        SRotation_single_particle(part, -sin_z, cos_z);
        XYShift_single_particle(part, -dx, -dy);

        // Update dead particles
        if (!is_alive){

            LocalParticle_set_state(part, XC_LOST_ON_ABSORBER);

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
