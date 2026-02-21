// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_IMPACTS_H
#define XCOLL_IMPACTS_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU

#include <xtrack/headers/track.h>


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
            InteractionRecordData_set__inter(record, i_slot, interaction);

            double charge_ratio = LocalParticle_get_charge_ratio(parent);
            double mass_ratio = charge_ratio / LocalParticle_get_chi(parent);
            double energy = ( LocalParticle_get_ptau(parent) + 1 / LocalParticle_get_beta0(parent)
                             ) * mass_ratio * LocalParticle_get_p0c(parent);
            // All fields have to be written, or the arrays will not have the same length
            // TODO: maybe this is not true, as we are setting by slot index? Don't the arrays come pre-initialised?
            InteractionRecordData_set_id_before(record, i_slot, LocalParticle_get_particle_id(parent));
            InteractionRecordData_set_s_before(record,  i_slot, LocalParticle_get_s(parent));
            InteractionRecordData_set_x_before(record,  i_slot, LocalParticle_get_x(parent));
            InteractionRecordData_set_px_before(record, i_slot, LocalParticle_get_px(parent));
            InteractionRecordData_set_y_before(record,  i_slot, LocalParticle_get_y(parent));
            InteractionRecordData_set_py_before(record, i_slot, LocalParticle_get_py(parent));
            InteractionRecordData_set_zeta_before(record,   i_slot, LocalParticle_get_zeta(parent));
            InteractionRecordData_set_delta_before(record,  i_slot, LocalParticle_get_delta(parent));
            InteractionRecordData_set_energy_before(record, i_slot, energy);
            InteractionRecordData_set_mass_before(record,   i_slot, mass_ratio*LocalParticle_get_mass0(parent));
            InteractionRecordData_set_charge_before(record, i_slot, charge_ratio*LocalParticle_get_q0(parent));
            // TODO: particle info
            InteractionRecordData_set_z_before(record, i_slot, -1);
            InteractionRecordData_set_a_before(record, i_slot, -1);
            InteractionRecordData_set_pdgid_before(record, i_slot, -1);

            // TODO: maybe this is not needed
            InteractionRecordData_set_id_after(record, i_slot, -1);
            InteractionRecordData_set_s_after(record, i_slot, -1);
            InteractionRecordData_set_x_after(record, i_slot, -1);
            InteractionRecordData_set_px_after(record, i_slot, -1);
            InteractionRecordData_set_y_after(record, i_slot, -1);
            InteractionRecordData_set_py_after(record, i_slot, -1);
            InteractionRecordData_set_zeta_after(record, i_slot, -1);
            InteractionRecordData_set_delta_after(record, i_slot, -1);
            InteractionRecordData_set_energy_after(record, i_slot, -1);
            InteractionRecordData_set_mass_after(record, i_slot, -1);
            InteractionRecordData_set_charge_after(record, i_slot, -1);
            InteractionRecordData_set_z_after(record, i_slot, -1);
            InteractionRecordData_set_a_after(record, i_slot, -1);
            InteractionRecordData_set_pdgid_after(record, i_slot, -1);
        }
    }
//     printf("Logging %i in slot %i\n", interaction, i_slot);
    return i_slot;
}

/*gpufun*/
void InteractionRecordData_log_child(InteractionRecordData record, int64_t i_slot, LocalParticle* child){
    if (record && i_slot>=0){
        double charge_ratio = LocalParticle_get_charge_ratio(child);
        double mass_ratio = charge_ratio / LocalParticle_get_chi(child);
        double energy = ( LocalParticle_get_ptau(child) + 1 / LocalParticle_get_beta0(child)
                         ) * mass_ratio * LocalParticle_get_p0c(child);
        InteractionRecordData_set_id_after(record, i_slot, LocalParticle_get_particle_id(child));
        InteractionRecordData_set_s_after(record,  i_slot, LocalParticle_get_s(child));
        InteractionRecordData_set_x_after(record,  i_slot, LocalParticle_get_x(child));
        InteractionRecordData_set_px_after(record, i_slot, LocalParticle_get_px(child));
        InteractionRecordData_set_y_after(record,  i_slot, LocalParticle_get_y(child));
        InteractionRecordData_set_py_after(record, i_slot, LocalParticle_get_py(child));
        InteractionRecordData_set_zeta_after(record,  i_slot, LocalParticle_get_zeta(child));
        InteractionRecordData_set_delta_after(record, i_slot, LocalParticle_get_delta(child));
        InteractionRecordData_set_energy_after(record, i_slot, energy);
        InteractionRecordData_set_mass_after(record,   i_slot, mass_ratio*LocalParticle_get_mass0(child));
        InteractionRecordData_set_charge_after(record, i_slot, charge_ratio*LocalParticle_get_q0(child));
        // TODO: particle info
        InteractionRecordData_set_z_after(record, i_slot, -1);
        InteractionRecordData_set_a_after(record, i_slot, -1);
        InteractionRecordData_set_pdgid_after(record, i_slot, -1);
//     printf("Slot %i: length %f\n", i_slot, ds);
    }
}

#endif /* XCOLL_IMPACTS_H */
