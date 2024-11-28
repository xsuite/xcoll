# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xpart as xp
import xpart.pdg as pdg
import xobjects as xo


def _drift(coll, particles, length):
    old_length = coll._equivalent_drift.length
    coll._equivalent_drift.length = length
    coll._equivalent_drift.track(particles)
    coll._equivalent_drift.length = old_length

def track(coll, particles):
    from ...beam_elements import FlukaCollimator
    if not isinstance(coll, FlukaCollimator):
        raise ValueError("Collimator is not a FlukaCollimator!\nCannot use FLUKA to track.")

    if not coll.active or not coll._tracking:
        _drift(coll, particles, coll.length)
        return

    npart = particles._num_active_particles
    if npart == 0:
        return

    # Check the server and whether it's initialised correctly
    from .engine import FlukaEngine
    engine = FlukaEngine().instance
    if not engine._flukaio_connected:
        raise ValueError(f"FlukaEngine not yet running!\nPlease do this first, by calling "
                       + f"xcoll.FlukaEngine.start(fluka_input_file.inp). "
                       + f"(id: {id(engine)})")

    if not engine._has_particle_ref:
        raise ValueError(f"FLUKA reference particle not set!\nPlease do this first, by calling "
                       + f"xcoll.FlukaEngine.set_particle_ref().")

    if 1.4*npart > engine._capacity:
        raise ValueError(f"Tracking {npart} particles but only {engine._capacity} allocated in "
                       + f"FlukaEngine!\nRemember to leave room for secondaries...")

    FlukaEngine.init_tracking(npart)

    if particles.particle_id.max() > engine.max_particle_id:
        raise ValueError(f"Some particles have an id that is higher than the highest id known "
                       + f"to FLUKA ({engine.max_particle_id}).\nThis could happen if this "
                       + f"particles object is larger than the first particles instance "
                       + f"tracked in this session, or if secondary particles are generated "
                       + f"somewhere else than FLUKA.\nIn that case, call "
                       + f"xcoll.FlukaEngine.init_tracking(max_particle_id) before tracking "
                       + f"with a value large enough to accommodate secondaries outside of FLUKA.\n"
                       + f"In any case, please stop and restart the FlukaEngine now.")

    if abs(particles.mass0 - engine.particle_ref.mass0) > 1e-3:
        raise ValueError("Error in reference mass of `particles`: not in sync with FLUKA reference particle!\n"
                       + "Rebuild the particles object using the FLUKA reference particle.")
    if abs(particles.q0 - engine.particle_ref.q0) > 1e-3:
        raise ValueError("Error in reference charge of `particles`: not in sync with FLUKA reference particle!\n"
                       + "Rebuild the particles object using the FLUKA reference particle.")
    if np.any([pdg_id == 0 for pdg_id in particles.pdg_id]):
        raise ValueError("Some particles are missing the pdg_id!")

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
    from .engine import FlukaEngine
    max_part = FlukaEngine.instance._capacity
    return np.concatenate((arr, np.zeros(max_part-arr.size, dtype=dtype)))


def track_core(coll, part):
    npart = part._num_active_particles
    from .engine import FlukaEngine
    engine = FlukaEngine.instance
    try:
        from .pyflukaf import track_fluka
    except ImportError as error:
        engine._warn_pyfluka(error)
        return

    max_part       = engine._capacity
    alive_at_entry = part.state > 0
    max_id         = part.particle_id[alive_at_entry].max()
    assert alive_at_entry.sum() == npart

    # Get particle data
    mass       = part.mass0*part.charge_ratio[alive_at_entry].copy() / part.chi[alive_at_entry].copy()
    charge     = part.q0*part.charge_ratio[alive_at_entry].copy()
    pdg_id     = part.pdg_id[alive_at_entry].copy()
    _, A, Z, _ = pdg.get_properties_from_pdg_id(pdg_id.copy())

    # Prepare arrays for FORTRAN
    # TODO: exact angles and drifts
    data = {}
    data['x']      = _expand(part.x[alive_at_entry].copy() * 1000.)
    data['xp']     = _expand(part.px[alive_at_entry].copy() * part.rpp[alive_at_entry].copy() * 1000.)
    data['y']      = _expand(part.y[alive_at_entry].copy() * 1000.)
    data['yp']     = _expand(part.py[alive_at_entry].copy() * part.rpp[alive_at_entry].copy() * 1000.)
    data['zeta']   = _expand(part.zeta[alive_at_entry].copy() * 1000.)
    data['e']      = _expand(part.energy[alive_at_entry].copy() / 1.e6)
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
    old_pid        = part.particle_id[alive_at_entry].copy()
    old_ppid       = part.parent_particle_id[alive_at_entry].copy()
    data['weight'] = _expand(part.weight[alive_at_entry].copy())
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
    if coll._acc_ionisation_loss < 0:
        coll._acc_ionisation_loss = 0.

    # send to fluka
    track_fluka(turn=turn_in+1,                # Turn indexing start from 1 with FLUKA IO (start from 0 with xpart)
                fluka_id=coll.fluka_id,
                length=coll.length + coll.length_front + coll.length_back,
                part_p0c=part.p0c[0],
                part_e0=part.energy0[0],
                alive_part=npart,
                max_part=max_part,
                x_part=data['x'],
                xp_part=data['xp'],
                y_part=data['y'],
                yp_part=data['yp'],
                zeta_part=data['zeta'],
                e_part=data['e'],
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
    new_ppid[mask_to_restore] = old_ppid[idx_to_restore]
    # When the parent changed (interaction), the parent ID must have been an alive particle at entry
    assert np.all([ppid in old_pid for ppid in new_ppid[~mask_to_restore]])

    # TODO: Impact Table
    #     Absorbed: trivial
    #     Not hit:    pid+ppid did not change and dE = 0
    #     Hit (MCS):  pid+ppid did not change and dE != 0
    #     Hit (nucl): above mask + children

    # IMPORTANT: we assume that a parent can never continue existing after an interaction. TODO: is this correct?

    # TODO: FLUKA returns particles that have undergone a nuclear interaction as new particles: we probably want to correct (after registering interaction)

    # Update existing particles  (these missed the collimator or only underwent elastic interactions)
    # ===============================================================================================
    mask_existing = new_pid <= max_id

    if np.any(mask_existing):
        idx_old  = np.array([np.where(part.particle_id==idx)[0][0] for idx in new_pid[mask_existing]])  # list of indices

        # Sanity check
        assert np.all(part.particle_id[idx_old] == new_pid[mask_existing])
        assert np.all(part.parent_particle_id[idx_old] == new_ppid[mask_existing])
        assert np.all(part.state[idx_old] > 0)

        # Update energy   TODO: in principle FLUKA uses pc, not energy!!
        E_diff = np.zeros(len(part.x))
        E_diff[idx_old] = part.energy[idx_old] - data['e'][:npart][mask_existing]*1.e6
        part.add_to_energy(-E_diff) # TODO: need to correct for weight
        coll._acc_ionisation_loss += np.sum(E_diff[idx_old])
        rpp = part.rpp[idx_old]    # This is now already updated by the new energy

        # Update other fields
        # TODO: exact angles etc
        part.x[idx_old]            = data['x'][:npart][mask_existing] / 1000.
        part.px[idx_old]           = data['xp'][:npart][mask_existing] / rpp / 1000.
        part.y[idx_old]            = data['y'][:npart][mask_existing] / 1000.
        part.py[idx_old]           = data['yp'][:npart][mask_existing] / rpp / 1000.
        part.zeta[idx_old]         = data['zeta'][:npart][mask_existing] / 1000.
        part.charge_ratio[idx_old] = data['q'][:npart][mask_existing] / part.q0
        part.chi[idx_old]          = data['q'][:npart][mask_existing] / (data['m'][:npart][mask_existing]*1.e6) \
                                               * part.mass0 / part.q0
        part.pdg_id[idx_old]       = data['pdg_id'][:npart][mask_existing]
        part.weight[idx_old]       = data['weight'][:npart][mask_existing]
        part.s[idx_old]            = s_in + coll.length + coll.length_front + coll.length_back

    # Little hack to set the dead particles, as idx_old is not a mask (but a list of indices)
    # (hence we cannot use ~idx_old)
    part.state[alive_at_entry]     = -334  # XC_LOST_ON_FLUKA
    if np.any(mask_existing):
        part.state[idx_old]        = 1     # These actually survived

    # Add new particles
    # ================
    mask_new = new_pid > max_id

    if not np.any(mask_new):
        part.reorganize()

    else:
        # Check that there is enough room in the particles object
        num_free = part._capacity - part._num_lost_particles - part._num_active_particles
        num_needed = np.sum(mask_new)
        if num_free < num_needed:
            raise RuntimeError(f"Too many particles generated by FLUKA ({num_needed} needed, "
                             + f"but only {num_free} free in particles object)!")
            # extra_capacity = int(1.2*(num_needed - num_free))  # 20% margin
            # # TODO: this does not work!!
            # part = xp.Particles.from_dict(part.to_dict(), _capacity=part._capacity+extra_capacity)

        # Sanity check
        idx_parents = np.array([np.where(part.particle_id==idx)[0][0] for idx in new_ppid[mask_new]])
        assert np.all(part.state[idx_parents] == -334)  # XC_LOST_ON_FLUKA

        # Create new particles
        # TODO: FLUKA uses pc; then calculate delta from p
        m     = data['m'][:npart][mask_new] * 1.e6
        E     = data['e'][:npart][mask_new] * 1.e6
        E0    = part.energy0[0]
        b0    = part.beta0[0]
        m0    = part.mass0
        delta = np.sqrt(1./(b0**2) * E**2/(E0**2) * m0**2/(m**2) - 1./(b0**2) + 1.) - 1.
        rpp   = 1. / (1. + delta)
        new_part = xp.Particles(_context=part._buffer.context,
                p0c = part.p0c[0],
                mass0 = part.mass0,
                q0 = part.q0,
                s = s_in + coll.length + coll.length_front + coll.length_back,
                x = data['x'][:npart][mask_new] / 1000.,
                px = data['xp'][:npart][mask_new] / rpp / 1000.,
                y = data['y'][:npart][mask_new] / 1000.,
                py = data['yp'][:npart][mask_new] / rpp / 1000.,
                zeta = data['zeta'][:npart][mask_new] / 1000.,
                delta = delta,
                mass_ratio = m / m0,
                charge_ratio = data['q'][:npart][mask_new] / part.q0,
                at_element = ele_in,
                at_turn = turn_in,
                pdg_id = data['pdg_id'][:npart][mask_new],
                parent_particle_id = new_ppid[mask_new],
                weight = data['weight'][:npart][mask_new],
                start_tracking_at_element = part.start_tracking_at_element)

        # Correct energy of parent particles: not everything was lost there
        # TODO: in principle FLUKA uses pc, not energy!!
        E_diff = np.bincount(idx_parents, weights=-new_part.energy, minlength=len(part.x))
        part.add_to_energy(E_diff)   # TODO: need to correct for weight

        # Add new particles
        new_part._init_random_number_generator()
        part.add_particles(new_part)
        engine.max_particle_id = part.particle_id.max()
