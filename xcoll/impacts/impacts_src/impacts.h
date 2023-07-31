// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #


#ifndef XCOLL_IMPACTS_H
#define XCOLL_IMPACTS_H

/*gpufun*/
void CollimatorImpacts_track_local_particle(CollimatorImpactsData el, LocalParticle* part0) {
    kill_all_particles(part0, XC_ERR_INVALID_TRACK);
}


/*gpufun*/
void CollimatorImpactsData_log(CollimatorImpactsData record, RecordIndex record_index, LocalParticle* parent,
                               LocalParticle* child, double parent_ds, double child_ds, int64_t interaction){
    // Record impact data
    if (record){
        // Get a slot in the record (this is thread safe)
        int64_t i_slot = RecordIndex_get_slot(record_index);
        // The returned slot id is negative if record is NULL or if record is full

        if (i_slot>=0){
            CollimatorImpactsData_set_at_element(record, i_slot, LocalParticle_get_at_element(parent));
            CollimatorImpactsData_set_at_turn(record, i_slot, LocalParticle_get_at_turn(parent));
            CollimatorImpactsData_set_s(record, i_slot, LocalParticle_get_s(parent));
            CollimatorImpactsData_set__inter(record, i_slot, interaction);

            // All fields have to be written, or the arrays will not have the same length
            CollimatorImpactsData_set_parent_id(record, i_slot, LocalParticle_get_particle_id(parent));
            CollimatorImpactsData_set_parent_ds(record, i_slot, parent_ds);
            CollimatorImpactsData_set_parent_x(record, i_slot, -1);
            CollimatorImpactsData_set_parent_px(record, i_slot, -1);
            CollimatorImpactsData_set_parent_y(record, i_slot, -1);
            CollimatorImpactsData_set_parent_py(record, i_slot, -1);
            CollimatorImpactsData_set_parent_zeta(record, i_slot, -1);
            CollimatorImpactsData_set_parent_delta(record, i_slot, -1);
            CollimatorImpactsData_set_parent_energy(record, i_slot, -1);
            // TODO: particle info
            CollimatorImpactsData_set_parent_mass(record, i_slot, -1);
            CollimatorImpactsData_set_parent_charge(record, i_slot, 1);
            CollimatorImpactsData_set_parent_z(record, i_slot, -1);
            CollimatorImpactsData_set_parent_a(record, i_slot, -1);
            CollimatorImpactsData_set_parent_pdgid(record, i_slot, -1);

            CollimatorImpactsData_set_child_id(record, i_slot, LocalParticle_get_particle_id(child));
            CollimatorImpactsData_set_child_ds(record, i_slot, child_ds);
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

#endif /* XCOLL_IMPACTS_H */