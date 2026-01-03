# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xtrack as xt
import xtrack.particles.pdg as pdg

from ...constants import (LOST_ON_FLUKA_COLL, MASSLESS_OR_NEUTRAL,
                          VIRTUAL_ENERGY, HIT_ON_FLUKA_COLL)


def track_pre(coll, particles):
    import xcoll as xc

    # Initialize ionisation loss accumulation variable
    if coll._acc_ionisation_loss < 0:
        coll._acc_ionisation_loss = 0.

    if not xc.fluka.engine.assert_ready_to_track_or_skip(
    coll, particles, _necessary_attributes=['fluka_id']):
        return False  # Stop tracking

    if not xc.fluka.engine._flukaio_connected:
        raise ValueError(f"FlukaEngine not yet running!\nPlease do this first, by calling "
                       + f"xcoll.FlukaEngine.start(fluka_input_file.inp). "
                       + f"(id: {id(xc.fluka.engine)})")

    npart = particles._num_active_particles
    if 1.4*npart > xc.fluka.engine.capacity:
        raise ValueError(f"Tracking {npart} particles but only {xc.fluka.engine.capacity} allocated in "
                       + f"FlukaEngine!\nRemember to leave room for secondaries...")
    if xc.fluka.engine.relative_capacity == 2 and xc.fluka.engine.particle_ref.pdg_id[0] != 2212:
        xc.fluka.engine._print("Warning: relative_capacity is set to 2. This is "
                             + "probably not enough for anything except protons.")

    xc.fluka.engine.init_tracking(npart+particles._num_lost_particles)

    if particles.particle_id.max() > xc.fluka.engine.max_particle_id:
        raise ValueError(f"Some particles have an id ({particles.particle_id.max()}) "
                       + f"that is higher than the highest id known to FLUKA ("
                       + f"{xc.fluka.engine.max_particle_id}).\nThis could happen if this "
                       + "particles object is larger than the first particles instance "
                       + "tracked in this session, or if secondary particles are generated "
                       + "somewhere else than FLUKA.\nIn that case, call "
                       + "xcoll.FlukaEngine.init_tracking(max_particle_id) before tracking "
                       + "with a value large enough to accommodate secondaries outside of "
                       + "FLUKA.\nIn any case, please stop and restart the FlukaEngine now.")

    coll._drift(particles, -coll.length_front)
    return True  # Continue tracking


def track_post(coll, particles):
    coll._drift(particles, -coll.length_back)
    alive_states = np.unique(particles.state[particles.state > 0])
    assert len(alive_states) <= 1, f"Unexpected alive particle states after tracking: {alive_states}"
    if len(alive_states) == 1:
        assert alive_states[0] == 1, f"Unexpected alive particle state after tracking: {alive_states[0]}"


def _expand(arr, available_capacity, dtype=float):
    return np.concatenate((arr, np.zeros(available_capacity, dtype=dtype)))


def track_core(coll, part):
    import xcoll as xc
    try:
        from pyflukaf import track_fluka
    except (ModuleNotFoundError, ImportError) as error:
        xc.fluka.engine._warn_pyfluka(error)
        return

    send_to_fluka  = part.state == HIT_ON_FLUKA_COLL
    npart          = send_to_fluka.sum()
    max_id         = part.particle_id[part.state > -9999].max()
    assert npart  <= part._num_active_particles
    if npart == 0:
        return

    # Get particle data
    m0         = part.mass0
    q0         = part.q0
    p0c        = part.p0c[0]
    E0         = part.energy0[0]
    beta0      = part.beta0[0]
    gamma0     = part.gamma0[0]
    mass       = part.mass[send_to_fluka]
    charge     = part.charge[send_to_fluka]
    pdg_id     = part.pdg_id[send_to_fluka]
    _, A, Z, _ = pdg.get_properties_from_pdg_id(pdg_id)
    A          = np.array([0 if pdgid < 0 else aa for aa, pdgid in zip(A, pdg_id)])  # FLUKA treats antiprotons as A = 0 = Z
    Z          = np.array([0 if pdgid < 0 else zz for zz, pdgid in zip(Z, pdg_id)])
    precision  = p0c * 1.e-6  # To avoid numerical issues like negative enervy. Ideally this should be 2.22e-15
    # TODO: VERY VERY BAD PRECISION

    # Decide how much extra capacity to send to FLUKA
    min_capacity = 50
    available_capacity = xc.fluka.engine.capacity - max_id - 1
    if available_capacity < min_capacity:
        raise ValueError("Not enough capacity in particles to accomodate secondaries.")
    available_capacity = min(int(np.ceil((xc.fluka.engine.relative_capacity - 1)*npart)),
                         available_capacity)
    available_capacity = max(available_capacity, min_capacity) # Some minimum value to be safe

    # Prepare arrays for FORTRAN
    data = {}
    data['x']      = _expand(part.x[send_to_fluka] * 1000., available_capacity)
    data['xp']     = _expand(part.px[send_to_fluka] * part.rpp[send_to_fluka] * 1000.,
                             available_capacity)
    data['y']      = _expand(part.y[send_to_fluka] * 1000., available_capacity)
    data['yp']     = _expand(part.py[send_to_fluka] * part.rpp[send_to_fluka] * 1000.,
                             available_capacity)
    data['zeta']   = _expand(part.zeta[send_to_fluka] * 1000., available_capacity)
    data['e']      = _expand(part.energy[send_to_fluka] / 1.e6, available_capacity)
    data['m']      = _expand(mass / 1.e6, available_capacity)
    data['q']      = _expand(charge.astype(np.int16), available_capacity, dtype=np.int16)
    data['A']      = _expand(A.astype(np.int32), available_capacity, dtype=np.int32)
    data['Z']      = _expand(Z.astype(np.int32), available_capacity, dtype=np.int32)
    data['pdg_id'] = _expand(pdg_id.astype(np.int32), available_capacity, dtype=np.int32)
    # FLUKA is 1-indexed
    data['pid']    = _expand(part.particle_id[send_to_fluka].astype(np.int32) + 1,
                             available_capacity, dtype=np.int32)
    # FLUKA does not use a parent ID, but a primary ID (hence not the direct parent but the first impact)
    # After one passage, there is no difference between parent ID and primary ID, but when a child gets
    # children in a second passage, we cannot trace them to the correct parent (only to the correct grand-
    # parent). To accommodate this, we do not send the parent ID to FLUKA but keep it in the particles
    # object.
    data['ppid']   = data['pid'].copy()
    data['weight'] = _expand(part.weight[send_to_fluka], available_capacity)
    # TODO: Hard-coded spin (currently not used)
    data['spin_x'] = _expand(np.zeros(npart), available_capacity)
    data['spin_y'] = _expand(np.zeros(npart), available_capacity)
    data['spin_z'] = _expand(np.zeros(npart), available_capacity)

    # Change npart to np.array to make it writable, store some initial data
    npart    = np.array(npart, dtype=np.int64)
    s_in     = part.s[send_to_fluka][0]
    ele_in   = part.at_element[send_to_fluka][0]
    turn_in  = part.at_turn[send_to_fluka][0]
    start    = part.start_tracking_at_element  # TODO: is this needed?

    # send to fluka
    track_fluka(turn=turn_in+1,     # Turn indexing start from 1 with FLUKA IO (start from 0 with xpart)
                fluka_id=coll.fluka_id,
                length=coll.length + coll.length_front + coll.length_back,
                alive_part=npart,
                max_part=npart + available_capacity,
                x_part=data['x'],
                xp_part=data['xp'], # FLUKA uses director cosine. This is exactly equal to px / (1+delta)
                y_part=data['y'],
                yp_part=data['yp'],
                zeta_part=data['zeta'],
                e_part=data['e'],   # FlukaIO pretends this is momentum, but the flukaserver source proves it is energy
                m_part=data['m'],
                q_part=data['q'],
                a_part=data['A'],
                z_part=data['Z'],
                pdg_id_part=data['pdg_id'],
                part_id=data['pid'],
                parent_id=data['ppid'],
                part_weight=data['weight'],
                spin_x_part=data['spin_x'],
                spin_y_part=data['spin_y'],
                spin_z_part=data['spin_z']
    )

    # Careful with all the masking!
    # Double-mask assignment does not work, e.g. part.state[mask1][mask2] = 1 will do nothing...

    new_pid  = data['pid'][:npart] - 1   # return to python 0-index
    new_ppid = data['ppid'][:npart] - 1  # return to python 0-index

    # TODO: Impact Table
    #     Absorbed: trivial
    #     Not hit:    pid+ppid did not change and dE = 0
    #     Hit (MCS):  pid+ppid did not change and dE != 0
    #     Hit (nucl): above mask + children

    # TODO: FLUKA returns particles that have undergone a nuclear interaction as new particles: we probably want to correct (after registering interaction)

    E_diff = np.zeros(len(part.x))


    # Update existing particles  (these missed the collimator or only underwent elastic interactions)
    # ===============================================================================================
    mask_alive = new_pid <= max_id
    idx_old = np.argsort(part.particle_id)
    pid_old_sorted = part.particle_id[idx_old]
    pos = np.searchsorted(pid_old_sorted, new_pid[mask_alive])
    idx_alive = idx_old[pos]

    # Kill particles that died just now
    # We cannot do part.state[~idx_alive] = LOST_ON_FLUKA_COLL (as idx_alive is not a mask)
    # So instead, we flag all as dead, and correct the alive ones later
    part.state[send_to_fluka] = -LOST_ON_FLUKA_COLL # Do not kill yet to avoid issues with energy updating

    if np.any(mask_alive):
        # Sanity check
        assert np.all(part.particle_id[idx_alive] == new_pid[mask_alive])
        assert np.all(part.state[idx_alive] == HIT_ON_FLUKA_COLL)

        # Update energy
        E_diff[idx_alive] = part.energy[idx_alive] - data['e'][:npart][mask_alive] * 1.e6
        if np.any(E_diff < -precision):
            raise ValueError(f"FLUKA returned particle with energy higher than incoming particle!")
        E_diff[E_diff < precision] = 0. # Lower cut on energy loss
        part.add_to_energy(-E_diff)
        part.weight[idx_alive]     = data['weight'][:npart][mask_alive]
        coll._acc_ionisation_loss += np.sum(E_diff[idx_alive]*part.weight[idx_alive])

        rpp = part.rpp[idx_alive]  # This is now already updated by the new energy
        part.x[idx_alive]          = data['x'][:npart][mask_alive] / 1000.
        part.px[idx_alive]         = data['xp'][:npart][mask_alive] / rpp / 1000. # This is exact because FLUKA uses director cosine
        part.y[idx_alive]          = data['y'][:npart][mask_alive] / 1000.
        part.py[idx_alive]         = data['yp'][:npart][mask_alive] / rpp / 1000. # This is exact because FLUKA uses director cosine
        part.zeta[idx_alive]       = data['zeta'][:npart][mask_alive] / 1000.
        part.state[idx_alive]      = 1  # These actually survived


    # Add new particles
    # ================
    mask_new = new_pid > max_id
    q_new = data['q'][:npart]
    pdg_id = data['pdg_id'][:npart]
    mask_new &= xc.fluka.engine._mask_particle_return_types(pdg_id, q_new)

    if np.any(mask_new):
        # Check that there is enough room in the particles object
        num_assigned = part._num_lost_particles + part._num_active_particles
        num_free = part._capacity - num_assigned
        num_needed = mask_new.sum()
        if num_free < num_needed:
            raise RuntimeError(f"Too many particles generated by FLUKA ({num_needed} needed, "
                             + f"but only {num_free} free in particles object)!")

        # Parent particle IDs
        parents = new_ppid[mask_new]
        pids = part.particle_id[send_to_fluka]
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

        # Create new particles
        E = data['e'][:npart][mask_new] * 1.e6
        m = data['m'][:npart][mask_new] * 1.e6
        q = data['q'][:npart][mask_new]
        if np.any(E < -precision):
            raise ValueError(f"FLUKA returned particles with negative energy!")
        if np.any(m < -precision):
            raise ValueError(f"FLUKA returned particles with negative mass!")
        # TODO: we set massless particles to have half the energy in their mass.
        # To be adapted when Xsuite can handle neutral and massless particles.
        mask_massless = np.abs(m) < 1.e-12
        mask_neutral = np.abs(q) < 1.e-12
        m[mask_massless | mask_neutral] = E[mask_massless | mask_neutral]/np.sqrt(2)  # Half the energy is put in the mass
        q[mask_massless | mask_neutral] = q0
        delta = np.sqrt(E/m*E/m - 1)/beta0/gamma0 - 1
        rpp   = 1. / (1. + delta)
        new_part = xt.Particles(_context=part._buffer.context,
                p0c = p0c,
                mass0 = m0,
                q0 = q0,
                s = s_in + coll.length + coll.length_front + coll.length_back,
                x = data['x'][:npart][mask_new] / 1000.,
                px = data['xp'][:npart][mask_new] * (1. + delta) / 1000.,
                y = data['y'][:npart][mask_new] / 1000.,
                py = data['yp'][:npart][mask_new] * (1. + delta) / 1000.,
                zeta = data['zeta'][:npart][mask_new] / 1000.,
                delta = delta,
                mass_ratio = m / m0,
                charge_ratio = q / q0,
                at_element = ele_in,
                at_turn = turn_in,
                pdg_id = data['pdg_id'][:npart][mask_new],
                particle_id = new_pid[mask_new],
                parent_particle_id = new_ppid[mask_new],
                weight = data['weight'][:npart][mask_new],
                start_tracking_at_element = start)

        # Correct the deposited energy of parent particles: not everything was lost there.
        E_children = np.bincount(idx_parents, weights=new_part.energy, minlength=part._capacity)
        if np.any(E_children < -precision):
            raise ValueError(f"FLUKA returned particles with summed energy higher than parent particle!")
        # If the parent survived, this should not be done but the energy should be subtracted
        # from the accumulated ionisation loss (as it is accounted for by the child)
        mask_parent_survived = (part.state==1) & (E_children > 0)
        if np.any(mask_parent_survived):
            if np.any(E_diff[mask_parent_survived] - E_children[mask_parent_survived]) < -precision:
                raise ValueError(f"FLUKA returned children with a surviving parent, however, there "
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
        new_part.state[mask_massless | mask_neutral] = -MASSLESS_OR_NEUTRAL
        part.add_particles(new_part)
        max_particle_id = new_pid.max()
        if max_particle_id <= xc.fluka.engine.max_particle_id:
            raise ValueError(f"FLUKA returned new particles with IDs {max_particle_id} that are "
                           + f"lower than the highest ID known ({xc.fluka.engine.max_particle_id}).\n"
                           + "This should not happen. Please report this issue to the developers.")
        xc.fluka.engine._max_particle_id = max_particle_id

    # Kill all flagged particles
    part.state[part.state==-LOST_ON_FLUKA_COLL] = LOST_ON_FLUKA_COLL
    part.state[part.state==-VIRTUAL_ENERGY] = VIRTUAL_ENERGY
    # TODO: we instantly kill massless or neutral particles as Xsuite is not ready to handle them.
    part.state[part.state==-MASSLESS_OR_NEUTRAL] = MASSLESS_OR_NEUTRAL

    # Reshuffle
    part.reorganize()

    # Give all particles the same s position to avoid numerical differences
    part.s[part.state > 0] = s_in + coll.length + coll.length_front + coll.length_back
