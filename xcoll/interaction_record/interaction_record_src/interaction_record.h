// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_IMPACTS_H
#define XCOLL_IMPACTS_H

#ifdef XO_CONTEXT_CPU
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU

#define XC_IR_SET_UNGUARDED(record, field, slot, value) \
    InteractionRecordData_set_##field(record, slot, value)

#define XC_IR_SET_GUARDED(record, field, slot, value) \
    if ((slot) < InteractionRecordData_len_##field(record)) { \
        InteractionRecordData_set_##field(record, slot, value); \
    }

#define XC_IR_LOG_PARENT_OPTIONAL(SET, record, slot, parent, energy, mass_ratio, pdgid) \
    SET(record, particle_id_before, slot, LocalParticle_get_particle_id(parent)); \
    SET(record, s_before, slot, LocalParticle_get_s(parent)); \
    SET(record, x_before, slot, LocalParticle_get_x(parent)); \
    SET(record, px_before, slot, LocalParticle_get_px(parent)); \
    SET(record, y_before, slot, LocalParticle_get_y(parent)); \
    SET(record, py_before, slot, LocalParticle_get_py(parent)); \
    SET(record, zeta_before, slot, LocalParticle_get_zeta(parent)); \
    SET(record, delta_before, slot, LocalParticle_get_delta(parent)); \
    SET(record, energy_before, slot, energy); \
    SET(record, mass_before, slot, mass_ratio*LocalParticle_get_mass0(parent)); \
    SET(record, pdgid_before, slot, pdgid); \
    SET(record, particle_id_after, slot, -1); \
    SET(record, s_after, slot, -1); \
    SET(record, x_after, slot, -1); \
    SET(record, px_after, slot, -1); \
    SET(record, y_after, slot, -1); \
    SET(record, py_after, slot, -1); \
    SET(record, zeta_after, slot, -1); \
    SET(record, delta_after, slot, -1); \
    SET(record, energy_after, slot, -1); \
    SET(record, mass_after, slot, -1); \
    SET(record, pdgid_after, slot, -1)

#define XC_IR_LOG_CHILD_OPTIONAL(SET, record, slot, child, energy, mass_ratio, pdgid) \
    SET(record, particle_id_after, slot, LocalParticle_get_particle_id(child)); \
    SET(record, s_after, slot, LocalParticle_get_s(child)); \
    SET(record, x_after, slot, LocalParticle_get_x(child)); \
    SET(record, px_after, slot, LocalParticle_get_px(child)); \
    SET(record, y_after, slot, LocalParticle_get_y(child)); \
    SET(record, py_after, slot, LocalParticle_get_py(child)); \
    SET(record, zeta_after, slot, LocalParticle_get_zeta(child)); \
    SET(record, delta_after, slot, LocalParticle_get_delta(child)); \
    SET(record, energy_after, slot, energy); \
    SET(record, mass_after, slot, mass_ratio*LocalParticle_get_mass0(child)); \
    SET(record, pdgid_after, slot, pdgid)


// TODO: do we need to pass RecordIndex?
// probably can do RecordIndex record_index = InteractionRecordData_getp__index(record);  ?
/*gpufun*/
int64_t InteractionRecordData_log(InteractionRecordData record, RecordIndex record_index, LocalParticle* parent,
                                  int64_t interaction, int64_t shape_id){
    // This can be used for a point-like interaction where there is no child (or because it's equal to the parent)
    // or to log the parent first, to be followed up with InteractionRecordData_log_child on the same slot

    int64_t i_slot = -1;
    if (record){
        // Get a slot in the record (this is thread safe)
        i_slot = RecordIndex_get_slot(record_index);
        // The returned slot id is negative if record is NULL or if record is full

        if (i_slot>=0){
            InteractionRecordData_set_at_turn(record, i_slot, LocalParticle_get_at_turn(parent));
            InteractionRecordData_set_at_element(record, i_slot, LocalParticle_get_at_element(parent));
            InteractionRecordData_set__inter(record, i_slot, interaction);
            InteractionRecordData_set_shape_id(record, i_slot, shape_id);

            double charge_ratio = LocalParticle_get_charge_ratio(parent);
            double mass_ratio = charge_ratio / LocalParticle_get_chi(parent);
            int64_t pdgid = LocalParticle_get_pdg_id(parent);
            double energy = ( LocalParticle_get_ptau(parent) + 1 / LocalParticle_get_beta0(parent)
                             ) * mass_ratio * LocalParticle_get_p0c(parent);
            if (InteractionRecordData_get__record_all_columns(record)) {
                XC_IR_LOG_PARENT_OPTIONAL(XC_IR_SET_UNGUARDED, record, i_slot,
                                          parent, energy, mass_ratio, pdgid);
            } else {
                XC_IR_LOG_PARENT_OPTIONAL(XC_IR_SET_GUARDED, record, i_slot,
                                          parent, energy, mass_ratio, pdgid);
            }
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
        int64_t pdgid = LocalParticle_get_pdg_id(child);
        if (InteractionRecordData_get__record_all_columns(record)) {
            XC_IR_LOG_CHILD_OPTIONAL(XC_IR_SET_UNGUARDED, record, i_slot,
                                     child, energy, mass_ratio, pdgid);
        } else {
            XC_IR_LOG_CHILD_OPTIONAL(XC_IR_SET_GUARDED, record, i_slot,
                                     child, energy, mass_ratio, pdgid);
        }
//     printf("Slot %i: length %f\n", i_slot, ds);
    }
}

#undef XC_IR_SET_UNGUARDED
#undef XC_IR_SET_GUARDED
#undef XC_IR_LOG_PARENT_OPTIONAL
#undef XC_IR_LOG_CHILD_OPTIONAL

#endif /* XCOLL_IMPACTS_H */
