# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt
import xtrack.particles.pdg as pdg

from ...headers.particle_states import LOST_ON_FLUKA_COLL, MASSLESS_OR_NEUTRAL, VIRTUAL_ENERGY


def _drift(coll, particles, length):
    old_length = coll._equivalent_drift.length
    coll._equivalent_drift.length = length
    coll._equivalent_drift.track(particles)
    coll._equivalent_drift.length = old_length

def track(coll, particles):
    import xcoll as xc

    # Initialize ionisation loss accumulation variable
    if coll._acc_ionisation_loss < 0:
        coll._acc_ionisation_loss = 0.

    if xc.fluka.engine.assert_ready_to_track_or_skip(coll, particles,
                                _necessary_attributes=['fluka_id']):
        return

    if not xc.fluka.engine._flukaio_connected:
        raise ValueError(f"FlukaEngine not yet running!\nPlease do this first, by calling "
                       + f"xcoll.FlukaEngine.start(fluka_input_file.inp). "
                       + f"(id: {id(xc.fluka.engine)})")

    npart = particles._num_active_particles
    if 1.4*npart > xc.fluka.engine.capacity:
        raise ValueError(f"Tracking {npart} particles but only {xc.fluka.engine.capacity} allocated in "
                       + f"FlukaEngine!\nRemember to leave room for secondaries...")

    xc.fluka.engine.init_tracking(npart)

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

    _drift(coll, particles, -coll.length_front)
    if coll.co is not None: # FLUKA collimators are centered; need to shift
        dx = coll.co[1][0]
        dy = coll.co[1][1]
        particles.x -= dx
        particles.y -= dy
    track_core(coll, particles)
    if coll.co is not None:
        particles.x += dx
        particles.y += dy
    _drift(coll, particles, -coll.length_back)


def _expand(arr, dtype=float):
    import xcoll as xc
    max_part = xc.fluka.engine.capacity
    return np.concatenate((arr, np.zeros(max_part-arr.size, dtype=dtype)))


def track_core(coll, part):
    import xcoll as xc
    npart = part._num_active_particles
    try:
        from pyflukaf import track_fluka
    except (ModuleNotFoundError, ImportError) as error:
        xc.fluka.engine._warn_pyfluka(error)
        return

    max_part       = xc.fluka.engine.capacity
    alive_at_entry = part.state > 0
    max_id         = part.particle_id[alive_at_entry].max()
    assert alive_at_entry.sum() == npart

    # Get particle data
    m0         = part.mass0
    q0         = part.q0
    p0c        = part.p0c[0]
    E0         = part.energy0[0]
    beta0      = part.beta0[0]
    gamma0     = part.gamma0[0]
    mass       = m0*part.charge_ratio[alive_at_entry] / part.chi[alive_at_entry]
    charge     = q0*part.charge_ratio[alive_at_entry]
    pdg_id     = part.pdg_id[alive_at_entry]
    _, A, Z, _ = pdg.get_properties_from_pdg_id(pdg_id)
    A          = np.array([0 if pdgid < 0 else aa for aa, pdgid in zip(A, pdg_id)])  # FLUKA treats antiprotons as A = 0 = Z
    Z          = np.array([0 if pdgid < 0 else zz for zz, pdgid in zip(Z, pdg_id)])
    precision  = p0c * 2.22e-15  # To avoid numerical issues like negative enervy

    # Prepare arrays for FORTRAN
    data = {}
    data['x']      = _expand(part.x[alive_at_entry] * 1000.)
    data['xp']     = _expand(part.px[alive_at_entry] * part.rpp[alive_at_entry] * 1000.)
    data['y']      = _expand(part.y[alive_at_entry] * 1000.)
    data['yp']     = _expand(part.py[alive_at_entry] * part.rpp[alive_at_entry] * 1000.)
    data['zeta']   = _expand(part.zeta[alive_at_entry] * 1000.)
    data['e']      = _expand(part.energy[alive_at_entry] / 1.e6)
    data['m']      = _expand(mass / 1.e6)
    data['q']      = _expand(charge.astype(np.int16), dtype=np.int16)
    data['A']      = _expand(A.astype(np.int32), dtype=np.int32)
    data['Z']      = _expand(Z.astype(np.int32), dtype=np.int32)
    data['pdg_id'] = _expand(pdg_id.astype(np.int32), dtype=np.int32)
    # FLUKA is 1-indexed
    data['pid']    = _expand(part.particle_id[alive_at_entry].astype(np.int32) + 1, dtype=np.int32)
    # FLUKA does not use a parent ID, but a primary ID (hence not the direct parent but the first impact)
    # After one passage, there is no difference between parent ID and primary ID, but when a child gets
    # children in a second passage, we cannot trace them to the correct parent (only to the correct grand-
    # parent). To accommodate this, we do not send the parent ID to FLUKA but store it manually here. Then
    # in FLUKA everything looks like a first passage, and we can restore the correct info later.
    data['ppid']   = data['pid'].copy()
    old_pid        = part.particle_id[alive_at_entry]
    old_ppid       = part.parent_particle_id[alive_at_entry]
    data['weight'] = _expand(part.weight[alive_at_entry])
    # TODO: Hard-coded spin (currently not used)
    data['spin_x'] = _expand(np.zeros(npart))
    data['spin_y'] = _expand(np.zeros(npart))
    data['spin_z'] = _expand(np.zeros(npart))

    # Change npart to np.array to make it writable, store some initial data
    npart    = np.array(npart, dtype=np.int64)
    s_in     = part.s[0]
    ele_in   = part.at_element[0]
    turn_in  = part.at_turn[0]
    start    = part.start_tracking_at_element  # TODO: is this needed?

    # send to fluka
    track_fluka(turn=turn_in+1,     # Turn indexing start from 1 with FLUKA IO (start from 0 with xpart)
                fluka_id=coll.fluka_id,
                length=coll.length + coll.length_front + coll.length_back,
                alive_part=npart,
                max_part=max_part,
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
    # Restore the parent IDs: if ppid != pid an interaction occurred and nothing needs to be done
    # (the real parent ID is the primary ID in this case). Otherwise, we restore the original parent ID
    mask_to_restore = new_pid == new_ppid
    idx_to_restore  = np.array([np.where(old_pid==idx)[0][0] for idx in new_pid[mask_to_restore]])
    if len(idx_to_restore) > 0:
        new_ppid[mask_to_restore] = old_ppid[idx_to_restore]
    # When the parent changed (interaction), the parent ID must have been an alive particle at entry
    assert np.all([ppid in old_pid for ppid in new_ppid[~mask_to_restore]])

    # TODO: Impact Table
    #     Absorbed: trivial
    #     Not hit:    pid+ppid did not change and dE = 0
    #     Hit (MCS):  pid+ppid did not change and dE != 0
    #     Hit (nucl): above mask + children

    # TODO: FLUKA returns particles that have undergone a nuclear interaction as new particles: we probably want to correct (after registering interaction)

    # Update existing particles  (these missed the collimator or only underwent elastic interactions)
    # ===============================================================================================
    mask_existing = new_pid <= max_id

    if np.any(mask_existing):
        # TODO: this is slooooow
        idx_old  = np.array([np.where(part.particle_id[alive_at_entry]==idx)[0][0]
                             for idx in new_pid[mask_existing]])  # list of indices

        # Sanity check
        assert np.all(part.particle_id[idx_old] == new_pid[mask_existing])
        assert np.all(part.parent_particle_id[idx_old] == new_ppid[mask_existing])
        assert np.all(part.state[idx_old] > 0)

        # Update momentum
        E_diff = np.zeros(len(part.x))
        E_diff[idx_old] = part.energy[idx_old] - data['e'][:npart][mask_existing] * 1.e6
        if np.any(E_diff < -precision):
            raise ValueError(f"FLUKA returned particle with energy higher than incoming particle!")
        E_diff[E_diff < precision] = 0. # Lower cut on energy loss
        part.add_to_energy(-E_diff)
        part.weight[idx_old]       = data['weight'][:npart][mask_existing]
        coll._acc_ionisation_loss += np.sum(E_diff[idx_old]*part.weight[idx_old])
        rpp = part.rpp[idx_old]    # This is now already updated by the new energy

        part.x[idx_old]            = data['x'][:npart][mask_existing] / 1000.
        part.px[idx_old]           = data['xp'][:npart][mask_existing] / rpp / 1000. # This is exact because FLUKA uses director cosine
        part.y[idx_old]            = data['y'][:npart][mask_existing] / 1000.
        part.py[idx_old]           = data['yp'][:npart][mask_existing] / rpp / 1000. # This is exact because FLUKA uses director cosine
        part.zeta[idx_old]         = data['zeta'][:npart][mask_existing] / 1000.
        part.charge_ratio[idx_old] = data['q'][:npart][mask_existing] / q0
        part.chi[idx_old]          = part.charge_ratio[idx_old] * m0 / (data['m'][:npart][mask_existing] * 1.e6)
        part.s[idx_old]            = s_in + coll.length + coll.length_front + coll.length_back
        part.pdg_id[idx_old]       = data['pdg_id'][:npart][mask_existing]

    # Little hack to set the dead particles, as idx_old is not a mask (but a list of indices)
    # (hence we cannot use ~idx_old)
    part.state[alive_at_entry]     = -LOST_ON_FLUKA_COLL # Do not kill yet to avoid issues with energy updating
    if np.any(mask_existing):
        part.state[idx_old]        = 1     # These actually survived

    # Add new particles
    # ================
    mask_new = new_pid > max_id

    if np.any(mask_new):
        # Check that there is enough room in the particles object
        num_assigned = part._num_lost_particles + part._num_active_particles
        num_free = part._capacity - num_assigned
        num_needed = mask_new.sum()
        if num_free < num_needed:
            raise RuntimeError(f"Too many particles generated by FLUKA ({num_needed} needed, "
                             + f"but only {num_free} free in particles object)!")

        idx_parents = np.array([np.where(part.particle_id[alive_at_entry]==idx)[0][0] for idx in new_ppid[mask_new]])

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
                px = data['xp'][:npart][mask_new] / rpp / 1000.,
                y = data['y'][:npart][mask_new] / 1000.,
                py = data['yp'][:npart][mask_new] / rpp / 1000.,
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
