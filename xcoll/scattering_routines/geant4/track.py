# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import io
import numpy as np

import xpart as xp
from xtrack.particles import LAST_INVALID_STATE

from ...constants import (LOST_WITHOUT_SPEC, LOST_ON_GEANT4_COLL,
    EXCITED_ION_STATE, MASSLESS_OR_NEUTRAL, VIRTUAL_ENERGY, HIT_ON_GEANT4_COLL)


# TODO: change collimasim to use xp, yp, energy/momentum, mass instead of px, py, delta, mass_ratio


# def _drift(coll, particles, length):
#     old_length = coll._equivalent_drift.length
#     coll._equivalent_drift.length = length
#     coll._equivalent_drift.track(particles)
#     coll._equivalent_drift.length = old_length

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
    send_to_geant4 = part.state < 1.e21 # Temporary: send all particles to Geant4
    # npart          = send_to_geant4.sum()
    # max_id         = part.particle_id[part.state > -9999].max()
    # assert npart  <= part._num_active_particles
    # if npart == 0:
    #     return

    # This temp delta is necessary because for primary particles, the coordinates are
    # modified in place. But for the longitudinal plane there are 3 coordinates that must
    # be updated, so pass a copy of the delta for the update in place and trigger the
    # correct update of the 3 coordinates later
    delta_temp = part._delta[send_to_geant4].copy()
    npart = part._num_active_particles + part._num_lost_particles

    q0 = part.q0
    m0 = part.mass0
    p0c = part.p0c[0]
    beta0 = part.beta0[0]
    s_in = part.s[send_to_geant4][0]
    precision  = p0c * 2.22e-15  # To avoid numerical issues like negative energy

    if xc.geant4.engine.reentry_protection_enabled:
        ### remove this part after geant4 bug fixed
        coords = {
            'x': part.x[send_to_geant4],
            'y': part.y[send_to_geant4],
            'px': part.px[send_to_geant4],
            'py': part.py[send_to_geant4],
            'zeta': part.zeta[send_to_geant4],
            'delta': delta_temp,
            'chi': part.chi[send_to_geant4],
            'charge_ratio': part.charge_ratio[send_to_geant4],
            's': part.s[send_to_geant4],
            'pdg_id': part.pdg_id[send_to_geant4],
            'particle_id': part.particle_id[send_to_geant4],
            'state': part.state[send_to_geant4],     # TODO:  does collimasim work if states are HIT_ON_GEANT4_COLL?
            'at_element': part.at_element[send_to_geant4],
            'at_turn': part.at_turn[send_to_geant4]
        }
        buf = io.BytesIO() # Use numpy.savez to serialize
        np.savez(buf, **coords)
        buf.seek(0)
        result_blob = xc.geant4.engine._g4link.add_particles_and_collimate_return(
                                                buf.getvalue(), coll.geant4_id)
        result_buf = io.BytesIO(result_blob) # Deserialize
        products = np.load(result_buf)
    else:
        coords = [part.x[send_to_geant4], part.y[send_to_geant4], part.px[send_to_geant4], part.py[send_to_geant4],
                  part.zeta[send_to_geant4], delta_temp, part.chi[send_to_geant4],
                  part.charge_ratio[send_to_geant4], part.s[send_to_geant4],
                  part.pdg_id[send_to_geant4], part.particle_id[send_to_geant4], part.state[send_to_geant4],
                  part.at_element[send_to_geant4], part.at_turn[send_to_geant4]]
        xc.geant4.engine._g4link.addParticles(coords)
        xc.geant4.engine._g4link.selectCollimator(coll.geant4_id)
        xc.geant4.engine._g4link.collimate()
        products = xc.geant4.engine._g4link.collimateReturn(coords)

    # Returned particles should not have changed type
    assert np.all(part.parent_particle_id[:npart] == products['parent_particle_id'][:npart])
    assert np.all(part.pdg_id[:npart] == products['pdg_id'][:npart])
    assert np.all(part.charge_ratio[:npart] == products['charge_ratio'][:npart])
    assert np.all(part.mass_ratio[:npart] == products['mass_ratio'][:npart])

    # Dead particles should not have changed anything
    was_dead = part.state[:npart] < 1
    assert np.all(part.x[:npart][was_dead]    == products['x'][:npart][was_dead])
    assert np.all(part.y[:npart][was_dead]    == products['y'][:npart][was_dead])
    assert np.all(part.px[:npart][was_dead]   == products['px'][:npart][was_dead])
    assert np.all(part.py[:npart][was_dead]   == products['py'][:npart][was_dead])
    assert np.all(part.zeta[:npart][was_dead] == products['zeta'][:npart][was_dead])
    assert np.all(part.s[:npart][was_dead]    == products['s'][:npart][was_dead])

    # Update particles that died just now
    returned_dead = (part.state[:npart] > 0) & (products['state'][:npart] == LOST_WITHOUT_SPEC)
    was_dead_but_untreated = (part.state[:npart] < 1) & (products['state'][:npart] == LOST_WITHOUT_SPEC)
    assert was_dead_but_untreated.sum() == 0
    part.at_turn[:npart][returned_dead]    = products['at_turn'][:npart][returned_dead]
    part.at_element[:npart][returned_dead] = products['at_element'][:npart][returned_dead]
    # TODO: check if different, and if yes, also store ionisation loss (can join mask with the below one)
    part.x[:npart][returned_dead]          = products['x'][:npart][returned_dead]
    part.y[:npart][returned_dead]          = products['y'][:npart][returned_dead]
    part.px[:npart][returned_dead]         = products['px'][:npart][returned_dead]     # TODO: BDSIM uses expanded xp, can we control that here to use exact?
    part.py[:npart][returned_dead]         = products['py'][:npart][returned_dead]
    part.zeta[:npart][returned_dead]       = products['zeta'][:npart][returned_dead]
    part.s[:npart][returned_dead]          = products['s'][:npart][returned_dead]

    # Update particles that are still alive
    returned_alive = products['state'][:npart] > 0
    # Energy needs special treatment
    new_delta = products['delta'][:npart][returned_alive]
    mass_ratio = part.mass_ratio[:npart][returned_alive]
    new_energy = np.sqrt((1+new_delta)**2 * p0c**2 * mass_ratio**2 + m0**2)
    E_diff = np.zeros(len(part.x))
    E_diff[:npart][returned_alive] = part.energy[:npart][returned_alive] - new_energy
    if np.any(E_diff < -precision):
        raise ValueError(f"Geant4 returned particle with energy higher than incoming particle!")
    E_diff[E_diff < precision] = 0. # Lower cut on energy loss
    part.add_to_energy(-E_diff)
    coll._acc_ionisation_loss += np.sum(E_diff)

    part.x[:npart][returned_alive] = products['x'][:npart][returned_alive]
    part.y[:npart][returned_alive] = products['y'][:npart][returned_alive]
    part.px[:npart][returned_alive] = products['px'][:npart][returned_alive]     # TODO: BDSIM uses expanded xp, can we control that here to use exact?
    part.py[:npart][returned_alive] = products['py'][:npart][returned_alive]
    part.zeta[:npart][returned_alive] = products['zeta'][:npart][returned_alive]
    part.s[:npart][returned_alive] = products['s'][:npart][returned_alive]
    part.at_element[:npart][returned_alive] = products['at_element'][:npart][returned_alive]
    part.at_turn[:npart][returned_alive] = products['at_turn'][:npart][returned_alive]

    if products['x'] is None or products['x'][npart] == -9999:
        # No new particles created in Geant4
        part.reorganize()

    else:
        # Add new particles created in Geant4
        mask_new = products['state'][npart:] > LAST_INVALID_STATE

        # Check that there is enough room in the particles object
        num_assigned = part._num_lost_particles + part._num_active_particles
        num_free = part._capacity - num_assigned
        num_needed = mask_new.sum()
        if num_free < num_needed:
            raise RuntimeError(f"Too many particles generated by Geant4 ({num_needed} needed, "
                             + f"but only {num_free} free in particles object)!")

        # Parent particle IDs
        parents = products['parent_particle_id'][npart:][mask_new]
        assert np.all(parents >= 0)
        assert parents.max() <= part.particle_id[:npart].max()
        idx_parents = np.array([np.where(part.particle_id[send_to_geant4]==idx)[0][0] for idx in parents])

        # Mass
        mass_ratio = products['mass_ratio'][npart:][mask_new]
        if np.any(mass_ratio < -precision):
            raise ValueError(f"Geant4 returned particle with negative mass!")
        massless = np.abs(mass_ratio) < 1.e-12

        # Charge
        charge_ratio = products['charge_ratio'][npart:][mask_new]
        neutral = np.abs(charge_ratio) < 1.e-12
        charge_ratio[massless | neutral] = q0

        # Energy
        delta = products['delta'][npart:][mask_new]
        # TODO: does NOT work for massless particles because NaN in collimasim
        # TODO: need virtual particles treatment
        # energy = np.sqrt((1+delta)**2 * p0c**2 * mass_ratio**2 + m0**2)
        # # TODO: we set massless particles to have half the energy in their mass.
        # # To be adapted when Xsuite can handle neutral and massless particles.
        # m[mask_massless | mask_neutral] = E[mask_massless | mask_neutral]/np.sqrt(2)  # Half the energy is put in the mass
        # q[mask_massless | mask_neutral] = q0
        # delta = np.sqrt(E/m*E/m - 1)/beta0/gamma0 - 1

        # new_delta[massless | neutral] = 
        # rpp   = 1. / (1. + delta)

        new_part = xp.Particles(_context=part._buffer.context,
                p0c = part.p0c[0],
                mass0 = part.mass0,
                q0 = part.q0,
                s = products['s'][npart:][mask_new],
                x = products['x'][npart:][mask_new],
                px = products['px'][npart:][mask_new],
                y = products['y'][npart:][mask_new],
                py = products['py'][npart:][mask_new],
                zeta = products['zeta'][npart:][mask_new],
                delta = delta,
                mass_ratio = mass_ratio,
                charge_ratio = charge_ratio,
                at_element = products['at_element'][npart:][mask_new],
                at_turn = products['at_turn'][npart:][mask_new],
                parent_particle_id = products['parent_particle_id'][npart:][mask_new],
                pdg_id = products['pdg_id'][npart:][mask_new])

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
        # TODO: we instantly kill massless or neutral particles as Xsuite is not ready to handle them.
        new_part.state[massless | neutral] = -MASSLESS_OR_NEUTRAL

        # Set the state of excited ions - not supported (will fail in BDSIM when resent)
        mask = (new_part.pdg_id > 999999999) & ((new_part.pdg_id % 10) != 0)
        new_part.state[mask] = -EXCITED_ION_STATE

        part.add_particles(new_part)

    # Kill all flagged particles
    part.state[part.state==LOST_WITHOUT_SPEC] = LOST_ON_GEANT4_COLL
    part.state[part.state==-VIRTUAL_ENERGY] = VIRTUAL_ENERGY
    part.state[part.state==-EXCITED_ION_STATE] = EXCITED_ION_STATE
    part.state[part.state==-MASSLESS_OR_NEUTRAL] = MASSLESS_OR_NEUTRAL

    # Reshuffle
    part.reorganize()

    # Give all particles the same s position to avoid numerical differences
    part.s[part.state > 0] = s_in + coll.length
