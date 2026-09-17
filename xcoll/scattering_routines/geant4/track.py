# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import io
import numpy as np

from xtrack import Particles

from ...constants import (LOST_ON_MATERIAL, LOST_ON_MATERIAL_SEC,
                          LOST_WITHOUT_SPEC, EXCITED_ION_STATE,
                          MASSLESS_OR_NEUTRAL, VIRTUAL_ENERGY,
                          VIRTUAL_ENERGY_SEC, HIT_ON_GEANT4,
                          HIT_ON_GEANT4_SEC, SECONDARY_PARTICLE)


def track_pre(coll, particles):
    import xcoll as xc

    if not xc.geant4.engine.assert_ready_to_track_or_skip(
    coll, particles, _necessary_attributes=['geant4_id']):
        return False  # Stop tracking

    # Correct for 250nm BDSIM margin
    # TODO: need to understand this margin better. What is BDSIM doing?
    # Does this work with tilts?
    coll._drift(particles, -coll.length_front)

    return True  # Continue tracking


def track_post(coll, particles):
    coll._drift(particles, -coll.length_back)

    # Ensure no leftover states
    assert np.sum(particles.state==LOST_WITHOUT_SPEC) == 0
    alive_states = np.unique(particles.state[particles.state > 0])
    alive_states = alive_states[(alive_states > SECONDARY_PARTICLE) & (alive_states <= 399)]
    if any(alive_states):
        raise ValueError(f"After tracking through Geant4, some particles still "
                         f"have temporary hit states: {alive_states}. This "
                         f"should not happen.")


def track_core(coll, part):
    import xcoll as xc
    xc.geant4.engine._g4link.clearData() # Clear the old data - bunch particles and hits

    send_to_geant4 = (part.state == HIT_ON_GEANT4) | (part.state == HIT_ON_GEANT4_SEC)
    was_secondary  = part.state == HIT_ON_GEANT4_SEC
    npart          = send_to_geant4.sum()
    assert npart  <= part._num_active_particles
    if npart == 0:
        return

    num_sent = send_to_geant4.sum()
    idx_sent = np.arange(len(part.x))[send_to_geant4]

    # Get particle data
    q0 = part.q0
    m0 = part.mass0
    p0c = part.p0c[0]
    beta0 = part.beta0[0]
    gamma0 = part.gamma0[0]
    s_in = part.s[send_to_geant4][0]
    ele_in = part.at_element[send_to_geant4][0]
    turn_in = part.at_turn[send_to_geant4][0]
    precision  = p0c * 1.e-12  # To avoid numerical issues like negative energy. Ideally this should be 2.22e-15

    rpp  = part.rpp[send_to_geant4]
    x    = part.x[send_to_geant4]
    xp   = part.px[send_to_geant4] * rpp   # Director cosine
    y    = part.y[send_to_geant4]
    yp   = part.py[send_to_geant4] * rpp   # Director cosine
    zeta = part.zeta[send_to_geant4]
    p    = p0c * (1 + part.delta[send_to_geant4]) * part.mass_ratio[send_to_geant4]
    q    = q0 * part.charge_ratio[send_to_geant4]
    weight = part.weight[send_to_geant4]
    pdgid  = part.pdg_id[send_to_geant4]
    pid    = part.particle_id[send_to_geant4]

    if xc.geant4.engine.reentry_protection_enabled:
        ### remove this part after geant4 bug fixed
        coords = {'x': x, 'xp': xp, 'y': y, 'yp': yp, 'zeta': zeta, 'p': p,
                  'q': q, 'weight': weight, 'pdgid': pdgid, 'id': pid}
        buf = io.BytesIO() # Use numpy.savez to serialize
        np.savez(buf, **coords)
        buf.seek(0)
        result_blob = xc.geant4.engine._g4link.add_particles_and_collimate_return(
                            buf.getvalue(), coll.geant4_id, num_sent)
        result_buf = io.BytesIO(result_blob) # Deserialize
        products = np.load(result_buf)
    else:
        coords = [x, xp, y, yp, zeta, p, q, weight, pdgid, pid]
        xc.geant4.engine._g4link.addParticles(coords)
        xc.geant4.engine._g4link.selectCollimator(coll.geant4_id)
        xc.geant4.engine._g4link.collimate()
        products = xc.geant4.engine._g4link.collimateReturn(num_sent)

    # Careful with all the masking!
    # Double-mask assignment does not work, e.g. part.state[mask1][mask2] = 1 will do nothing...

    # Kill particles that died just now
    returned_dead = products['state'][:num_sent] == LOST_WITHOUT_SPEC
    idx_dead = idx_sent[returned_dead]
    part.state[idx_dead] = -LOST_ON_MATERIAL
    part.state[(part.state == -LOST_ON_MATERIAL) & was_secondary] = -LOST_ON_MATERIAL_SEC

    # Update particles that are still alive
    returned_alive = products['state'][:num_sent] > 0
    idx_alive = idx_sent[returned_alive]    # BDSIM keeps particle order
    idx_returned_alive = np.nonzero(returned_alive)[0]
    idx_alive_prim = idx_alive[part.state[idx_alive] == HIT_ON_GEANT4]
    idx_alive_sec  = idx_alive[part.state[idx_alive] == HIT_ON_GEANT4_SEC]

    # Energy needs special treatment
    m_in = part.mass[idx_alive]
    new_p = products['p'][idx_returned_alive]
    new_energy = np.sqrt(new_p**2 + m_in**2)
    E_diff = np.zeros(len(part.x))
    E_diff[idx_alive] = part.energy[idx_alive] - new_energy
    if np.any(E_diff < -precision):
        raise ValueError(f"Geant4 returned particle with energy higher than incoming particle!")
    E_diff[E_diff < precision] = 0. # Lower cut on energy loss
    part.add_to_energy(-E_diff)
    part.weight[idx_alive] = products['weight'][idx_returned_alive]
    coll._acc_ionisation_loss     += np.sum(E_diff[idx_alive_prim]*part.weight[idx_alive_prim])
    coll._acc_ionisation_loss_sec += np.sum(E_diff[idx_alive_sec]*part.weight[idx_alive_sec])
    mask_hit = ~np.isclose(E_diff, 0.) # Hit (MCS)
    # First mark the state of surviving particles as they were before
    part.state[idx_alive_prim] = 1
    part.state[idx_alive_sec] = SECONDARY_PARTICLE
    scattered_ids = part.particle_id[mask_hit]

    rpp  = part.rpp[idx_alive]
    part.x[idx_alive]    = products['x'][idx_returned_alive]
    part.px[idx_alive]   = products['xp'][idx_returned_alive] / rpp   # Director cosine back to px
    part.y[idx_alive]    = products['y'][idx_returned_alive]
    part.py[idx_alive]   = products['yp'][idx_returned_alive] / rpp   # Director cosine back to py
    part.zeta[idx_alive] = products['zeta'][idx_returned_alive]

    # Add new particles created in Geant4
    q_new = products['q'][num_sent:]
    pdg_id = products['pdg_id'][num_sent:]
    mask_new = xc.geant4.engine._mask_particle_return_types(pdg_id, q_new)
    idx_new = np.nonzero(mask_new)[0] + num_sent
    assert np.all(products['state'][idx_new] == 1)

    if np.any(mask_new):
        # Check that there is enough room in the particles object
        num_assigned = part._num_lost_particles + part._num_active_particles
        num_free = part._capacity - num_assigned
        num_needed = mask_new.sum()
        if num_free < num_needed:
            raise RuntimeError(f"Too many particles generated by Geant4 ({num_needed} needed, "
                             + f"but only {num_free} free in particles object)!")

        # Parent particle IDs
        parents = products['parent_particle_id'][idx_new]
        pids = part.particle_id[send_to_geant4]
        max_id_in = int(pids.max())
        assert np.all(parents >= 0)        # Sanity check
        assert parents.max() <= max_id_in  # Sanity check

        # Build dense lookup table: ID -> index in `pids`
        lookup = np.full(max_id_in + 1, -1, dtype=np.int64)
        lookup[pids] = np.arange(pids.size, dtype=np.int64)
        idx_parents_sent = lookup[parents]   # MUCH faster than solution with np.where
        if np.any(idx_parents_sent < 0):
            missing = np.unique(parents[idx_parents_sent < 0])
            raise RuntimeError(f"Parent IDs not found in particle_id: {missing}")
        idx_parents = idx_sent[idx_parents_sent]

        # Mass
        m_new = products['m'][idx_new]
        if np.any(m_new < -precision):
            raise ValueError(f"Geant4 returned particle with negative mass!")
        massless = np.abs(m_new) < 1.e-12

        # Charge
        q_new = products['q'][idx_new]
        neutral = np.abs(q_new) < 1.e-12
        massless_or_neutral = massless | neutral

        # Energy
        p_new = products['p'][idx_new]
        E_new = np.sqrt(p_new**2 + m_new**2)
        m_new[massless_or_neutral] = E_new[massless_or_neutral]/np.sqrt(2)
        q_new[massless_or_neutral] = q0
        delta = np.sqrt(E_new/m_new*E_new/m_new - 1)/beta0/gamma0 - 1

        new_state = np.ones(len(idx_new), dtype=np.int64)
        if coll.mark_scattered_particles:
            new_state[:] = SECONDARY_PARTICLE
        else:
            sec_states = np.array([-LOST_ON_MATERIAL_SEC, -VIRTUAL_ENERGY_SEC,
                                   SECONDARY_PARTICLE])
            parent_states = part.state[idx_parents]
            new_state[np.isin(parent_states, sec_states)] = SECONDARY_PARTICLE

        new_part = Particles(_context=part._buffer.context,
                p0c = part.p0c[0],
                mass0 = part.mass0,
                q0 = part.q0,
                s = s_in + coll.length,
                x = products['x'][idx_new],
                px = products['xp'][idx_new] * (1 + delta), # Director cosine back to px
                y = products['y'][idx_new],
                py = products['yp'][idx_new] * (1 + delta), # Director cosine back to py
                zeta = products['zeta'][idx_new],
                delta = delta,
                mass_ratio = m_new/m0,
                charge_ratio = q_new/q0,
                at_element = ele_in,
                at_turn = turn_in,
                state = new_state,
                parent_particle_id = products['parent_particle_id'][idx_new],
                pdg_id = products['pdg_id'][idx_new],
                weight = products['weight'][idx_new]
        )

        # Correct the deposited energy of parent particles: not everything was lost there.
        E_children = np.bincount(idx_parents, weights=new_part.energy, minlength=part._capacity)
        if np.any(E_children < -precision):
            raise ValueError(f"Geant4 returned particles with summed energy higher than parent particle!")
        # If the parent survived, this should not be done but the energy should be subtracted
        # from the accumulated ionisation loss (as it is accounted for by the child)
        mask_parent_survived_prim = (part.state==1) & (E_children > 0)
        mask_parent_survived_sec = (part.state==SECONDARY_PARTICLE) & (E_children > 0)
        mask_parent_survived = mask_parent_survived_prim | mask_parent_survived_sec
        if np.any(mask_parent_survived):
            if np.any(E_diff[mask_parent_survived] - E_children[mask_parent_survived] < -precision):
                raise ValueError(f"Geant4 returned children with a surviving parent, however, there "
                               + f"was a larger energy loss than the children energy!")
            coll._acc_ionisation_loss -= np.sum(E_children[mask_parent_survived_prim]
                                                * part.weight[mask_parent_survived_prim])
            coll._acc_ionisation_loss_sec -= np.sum(E_children[mask_parent_survived_sec]
                                                    * part.weight[mask_parent_survived_sec])
            E_children[mask_parent_survived] = 0.
        # If the deposited energy is lower than the rest mass, it cannot be represented by the original
        # particle. We make a virtual particle with fake mass (E=2m) to avoid negative square roots.
        mask_virtual = (part.energy - E_children < part.mass) & (E_children > 0)
        if np.any(mask_virtual):
            # import ipdb; ipdb.set_trace()
            virtual_mass = (part.energy - E_children)[mask_virtual] / np.sqrt(2)
            new_ptau = part.ptau.copy()
            old_mass = part.mass[mask_virtual]
            old_ptau = new_ptau[mask_virtual]
            part.chi[mask_virtual] = m0 / virtual_mass * part.charge_ratio[mask_virtual]
            new_ptau[mask_virtual] = (1/beta0 + old_ptau)*old_mass/virtual_mass - 1/beta0
            part.update_ptau(new_ptau)
            part.state[mask_virtual & (part.state==1)] = -VIRTUAL_ENERGY
            part.state[mask_virtual & (part.state==-LOST_ON_MATERIAL)] = -VIRTUAL_ENERGY
            part.state[mask_virtual & (part.state==-LOST_ON_MATERIAL_SEC)] = -VIRTUAL_ENERGY_SEC
            part.state[mask_virtual & (part.state==SECONDARY_PARTICLE)] = -VIRTUAL_ENERGY_SEC
        # Now update the parent energies
        part.add_to_energy(-E_children)

        # Add new particles
        new_part._init_random_number_generator()
        # TODO: we kill massless or neutral particles as Xsuite is not ready to handle them.
        new_part.state[massless_or_neutral] = -MASSLESS_OR_NEUTRAL

        # Set the state of excited ions - not supported (will fail in BDSIM when resent)
        # TODO: excited ion states are not supported but created by BDSIM
        mask = (new_part.pdg_id > 999999999) & ((new_part.pdg_id % 10) != 0)
        new_part.state[mask] = -EXCITED_ION_STATE

        part.add_particles(new_part)

    # Kill all flagged particles
    part.state[part.state==-LOST_ON_MATERIAL] = LOST_ON_MATERIAL
    part.state[part.state==-LOST_ON_MATERIAL_SEC] = LOST_ON_MATERIAL_SEC
    part.state[part.state==-VIRTUAL_ENERGY] = VIRTUAL_ENERGY
    part.state[part.state==-VIRTUAL_ENERGY_SEC] = VIRTUAL_ENERGY_SEC
    part.state[part.state==-EXCITED_ION_STATE] = EXCITED_ION_STATE
    # TODO: we instantly kill massless or neutral particles as Xsuite is not ready to handle them.
    part.state[part.state==-MASSLESS_OR_NEUTRAL] = MASSLESS_OR_NEUTRAL

    # Finally, mark scattered particles that survived, if requested.
    # This can only be done now to avoid issues with the energy updating
    # of virtual particles.
    if coll.mark_scattered_particles and len(scattered_ids) > 0:
        idx_sort = np.argsort(part.particle_id)
        pid_sorted = part.particle_id[idx_sort]
        pos = np.searchsorted(pid_sorted, scattered_ids)
        idx_scattered = idx_sort[pos]
        assert np.all(part.particle_id[idx_scattered] == scattered_ids)
        part.state[idx_scattered] = SECONDARY_PARTICLE

    # Reshuffle
    part.reorganize()

    # Give all particles the same s position to avoid numerical differences
    part.s[part.state > 0] = s_in + coll.length
