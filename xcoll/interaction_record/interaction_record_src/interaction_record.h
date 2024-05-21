// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #


#ifndef XCOLL_IMPACTS_H
#define XCOLL_IMPACTS_H

// TODO: do we need to pass RecordIndex?
// probably can do RecordIndex record_index = InteractionRecordData_getp__index(record);  ?
/*gpufun*/
int64_t InteractionRecordData_log(InteractionRecordData record, RecordIndex record_index, LocalParticle* parent,
                                  int64_t interaction){
    // This can be used for a point-like interaction where there is no child (or because it's equal to the parent)
    // or to log the parent first, to be followed up with InteractionRecordData_log_child on the same slot

    int64_t i_slot = -1;
    if (record){
        // Get a slot in the record (this is thread safe)
        i_slot = RecordIndex_get_slot(record_index);
        // The returned slot id is negative if record is NULL or if record is full

        if (i_slot>=0){
            InteractionRecordData_set_at_element(record, i_slot, LocalParticle_get_at_element(parent));
            InteractionRecordData_set_at_turn(record, i_slot, LocalParticle_get_at_turn(parent));
            InteractionRecordData_set_ds(record, i_slot, 0);
            InteractionRecordData_set__inter(record, i_slot, interaction);

            double charge_ratio = LocalParticle_get_charge_ratio(parent);
            double mass_ratio = charge_ratio / LocalParticle_get_chi(parent);
            double energy = ( LocalParticle_get_ptau(parent) + 1 / LocalParticle_get_beta0(parent)
                             ) * mass_ratio * LocalParticle_get_p0c(parent);
            // All fields have to be written, or the arrays will not have the same length
            // TODO: maybe this is not true, as we are setting by slot index? Don't the arrays come pre-initialised?
            InteractionRecordData_set_parent_id(record, i_slot, LocalParticle_get_particle_id(parent));
            InteractionRecordData_set_parent_x(record,  i_slot, LocalParticle_get_x(parent));
            InteractionRecordData_set_parent_px(record, i_slot, LocalParticle_get_px(parent));
            InteractionRecordData_set_parent_y(record,  i_slot, LocalParticle_get_y(parent));
            InteractionRecordData_set_parent_py(record, i_slot, LocalParticle_get_py(parent));
            InteractionRecordData_set_parent_zeta(record,   i_slot, LocalParticle_get_zeta(parent));
            InteractionRecordData_set_parent_delta(record,  i_slot, LocalParticle_get_delta(parent));
            InteractionRecordData_set_parent_energy(record, i_slot, energy);
            InteractionRecordData_set_parent_mass(record,   i_slot, mass_ratio*LocalParticle_get_mass0(parent));
            InteractionRecordData_set_parent_charge(record, i_slot, charge_ratio*LocalParticle_get_q0(parent));
            // TODO: particle info
            InteractionRecordData_set_parent_z(record, i_slot, -1);
            InteractionRecordData_set_parent_a(record, i_slot, -1);
            InteractionRecordData_set_parent_pdgid(record, i_slot, -1);

            // TODO: maybe this is not needed
            InteractionRecordData_set_child_id(record, i_slot, -1);
            InteractionRecordData_set_child_x(record, i_slot, -1);
            InteractionRecordData_set_child_px(record, i_slot, -1);
            InteractionRecordData_set_child_y(record, i_slot, -1);
            InteractionRecordData_set_child_py(record, i_slot, -1);
            InteractionRecordData_set_child_zeta(record, i_slot, -1);
            InteractionRecordData_set_child_delta(record, i_slot, -1);
            InteractionRecordData_set_child_energy(record, i_slot, -1);
            InteractionRecordData_set_child_mass(record, i_slot, -1);
            InteractionRecordData_set_child_charge(record, i_slot, -1);
            InteractionRecordData_set_child_z(record, i_slot, -1);
            InteractionRecordData_set_child_a(record, i_slot, -1);
            InteractionRecordData_set_child_pdgid(record, i_slot, -1);
        }
    }
//     printf("Logging %i in slot %i\n", interaction, i_slot);
    return i_slot;
}

/*gpufun*/
void InteractionRecordData_log_child(InteractionRecordData record, int64_t i_slot, LocalParticle* child, double ds){
    if (record && i_slot>=0){
        InteractionRecordData_set_ds(record, i_slot, ds);

        double charge_ratio = LocalParticle_get_charge_ratio(child);
        double mass_ratio = charge_ratio / LocalParticle_get_chi(child);
        double energy = ( LocalParticle_get_ptau(child) + 1 / LocalParticle_get_beta0(child)
                         ) * mass_ratio * LocalParticle_get_p0c(child);
        InteractionRecordData_set_child_id(record, i_slot, LocalParticle_get_particle_id(child));
        InteractionRecordData_set_child_x(record,  i_slot, LocalParticle_get_x(child));
        InteractionRecordData_set_child_px(record, i_slot, LocalParticle_get_px(child));
        InteractionRecordData_set_child_y(record,  i_slot, LocalParticle_get_y(child));
        InteractionRecordData_set_child_py(record, i_slot, LocalParticle_get_py(child));
        InteractionRecordData_set_child_zeta(record,  i_slot, LocalParticle_get_zeta(child));
        InteractionRecordData_set_child_delta(record, i_slot, LocalParticle_get_delta(child));
        InteractionRecordData_set_child_energy(record, i_slot, energy);
        InteractionRecordData_set_child_mass(record,   i_slot, mass_ratio*LocalParticle_get_mass0(child));
        InteractionRecordData_set_child_charge(record, i_slot, charge_ratio*LocalParticle_get_q0(child));
        // TODO: particle info
        InteractionRecordData_set_child_z(record, i_slot, -1);
        InteractionRecordData_set_child_a(record, i_slot, -1);
        InteractionRecordData_set_child_pdgid(record, i_slot, -1);
//     printf("Slot %i: length %f\n", i_slot, ds);
    }
}

#endif /* XCOLL_IMPACTS_H */
