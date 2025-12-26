# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import io
import numpy as np

from xtrack import Particles

from ...constants import (LOST_WITHOUT_SPEC, LOST_ON_GEANT4_COLL,
    EXCITED_ION_STATE, MASSLESS_OR_NEUTRAL, VIRTUAL_ENERGY, HIT_ON_GEANT4_COLL)


def track_pre(coll, particles):
    import xcoll as xc

    # Initialize ionisation loss accumulation variable
    if coll._acc_ionisation_loss < 0:
        coll._acc_ionisation_loss = 0.

    if not xc.geant4.engine.assert_ready_to_track_or_skip(
    coll, particles, _necessary_attributes=['geant4_id']):
        return False  # Stop tracking

    return True  # Continue tracking


def track_post(coll, particles):
    alive_states = np.unique(particles.state[particles.state > 0])
    assert len(alive_states) <= 1, f"Unexpected alive particle states after tracking: {alive_states}"
    if len(alive_states) == 1:
        assert alive_states[0] == 1, f"Unexpected alive particle state after tracking: {alive_states[0]}"


# TODO: need to rework this logic with HIT_ON_GEANT4_COLL


def track_core(coll, part):
    import xcoll as xc
    xc.geant4.engine._g4link.clearData() # Clear the old data - bunch particles and hits

    # send_to_geant4 = part.state == HIT_ON_GEANT4_COLL
    send_to_geant4 = part.state > 0
    # npart          = send_to_geant4.sum()
    # max_id         = part.particle_id[part.state > -9999].max()
    # assert npart  <= part._num_active_particles
    # if npart == 0:
    #     return

    num_sent = send_to_geant4.sum()

    coll._drift(part, -250e-9) # Margin before Geant4 collimator; added by BDSIM geometry

    # Get particle data
    q0 = part.q0
    m0 = part.mass0
    p0c = part.p0c[0]
    beta0 = part.beta0[0]
    s_in = part.s[send_to_geant4][0]
    ele_in = part.at_element[send_to_geant4][0]
    turn_in = part.at_turn[send_to_geant4][0]
    precision  = p0c * 1.e-12  # To avoid numerical issues like negative enervy. Ideally this should be 2.22e-15

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
    idx_dead = np.arange(len(part.x))[send_to_geant4][returned_dead]
    part.state[idx_dead] = -LOST_ON_GEANT4_COLL

    # Update particles that are still alive
    returned_alive = products['state'][:num_sent] > 0
    idx_alive = np.arange(len(part.x))[send_to_geant4][returned_alive]    # BDSIM keeps particle order
    idx_returned_alive = np.nonzero(returned_alive)[0]

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
    coll._acc_ionisation_loss += np.sum(E_diff[idx_alive]*part.weight[idx_alive])

    rpp  = part.rpp[idx_alive]
    part.x[idx_alive]      = products['x'][idx_returned_alive]
    part.px[idx_alive]     = products['xp'][idx_returned_alive] / rpp   # Director cosine back to px
    part.y[idx_alive]      = products['y'][idx_returned_alive]
    part.py[idx_alive]     = products['yp'][idx_returned_alive] / rpp   # Director cosine back to py
    part.zeta[idx_alive]   = products['zeta'][idx_returned_alive]

    # Add new particles created in Geant4
    q_new = products['q'][num_sent:]
    pdg_id = products['pdg_id'][num_sent:]
    mask_new = xc.geant4.engine._mask_particle_return_types(pdg_id, q_new)
    idx_new = np.nonzero(mask_new)[0] + num_sent
    assert np.all(products['state'][idx_new] == 1)

    if not np.any(mask_new):
        # No new particles created in Geant4
        part.reorganize()

    else:
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
        idx_parents = lookup[parents]   # MUCH faster than solution with np.where
        if np.any(idx_parents < 0):
            missing = np.unique(parents[idx_parents < 0])
            raise RuntimeError(f"Parent IDs not found in particle_id: {missing}")

        # Mass
        m_new = products['m'][idx_new]
        if np.any(m_new < -precision):
            raise ValueError(f"Geant4 returned particle with negative mass!")
        massless = np.abs(m_new) < 1.e-12
        m_new[massless] = m0  # TODO

        # Charge
        q_new = products['q'][idx_new]
        neutral = np.abs(q_new) < 1.e-12
        q_new[massless | neutral] = q0  # TODO

        # Energy
        p_new = products['p'][idx_new]
        delta = np.zeros_like(p_new)
        delta[massless] = p_new[massless]/p0c - 1
        delta[~massless] = p_new[~massless]/p0c * m0/m_new[~massless] - 1

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
        mask_parent_survived = (part.state==1) & (E_children > 0)
        if np.any(mask_parent_survived):
            if np.any(E_diff[mask_parent_survived] - E_children[mask_parent_survived]) < -precision:
                raise ValueError(f"Geant4 returned children with a surviving parent, however, there "
                               + f"was a larger energy loss than the children energy!")
            coll._acc_ionisation_loss -= np.sum(E_children[mask_parent_survived]*part.weight[mask_parent_survived])
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
            part.state[mask_virtual] = -VIRTUAL_ENERGY
        # Now update the parent energies
        part.add_to_energy(-E_children)

        # Add new particles
        new_part._init_random_number_generator()
        # TODO: we kill massless or neutral particles as Xsuite is not ready to handle them.
        new_part.state[massless | neutral] = -MASSLESS_OR_NEUTRAL

        # Set the state of excited ions - not supported (will fail in BDSIM when resent)
        mask = (new_part.pdg_id > 999999999) & ((new_part.pdg_id % 10) != 0)
        new_part.state[mask] = -EXCITED_ION_STATE

        part.add_particles(new_part)

    # Kill all flagged particles
    part.state[part.state==-LOST_ON_GEANT4_COLL] = LOST_ON_GEANT4_COLL
    part.state[part.state==-VIRTUAL_ENERGY] = VIRTUAL_ENERGY
    part.state[part.state==-EXCITED_ION_STATE] = EXCITED_ION_STATE
    part.state[part.state==-MASSLESS_OR_NEUTRAL] = MASSLESS_OR_NEUTRAL

    # Ensure no leftover states
    assert np.sum(part.state==LOST_WITHOUT_SPEC) == 0

    # Reshuffle
    part.reorganize()

    # Give all particles the same s position to avoid numerical differences
    part.s[part.state > 0] = s_in + coll.length
