# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xpart as xp
import xpart.pdg as pdg
import xobjects as xo


def drift_6d(particles, length):
    npart = particles._num_active_particles
    rpp = particles.rpp[:npart]
    xp = particles.px[:npart] * rpp
    yp = particles.py[:npart] * rpp
    dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
    particles.x[:npart] += xp * length
    particles.y[:npart] += yp * length
    particles.s[:npart] += length
    particles.zeta[:npart] += dzeta*length
    return


def track(collimator, particles):  # TODO: write impacts
    from ...beam_elements import FlukaCollimator
    if not isinstance(collimator, FlukaCollimator):
        raise ValueError("Collimator is not a FlukaCollimator!\nCannot use fluka to track.")

    if not collimator.active or not collimator._tracking:
        # Drift full length
        # TODO: when in C, drifting should call routine from drift element
        #       such that can be exact etc
        drift_6d(particles, collimator.length)
        return

    npart = particles._num_active_particles
    if npart == 0:
        return

    # Check the server and whether it's initialised correctly
    from .engine import FlukaEngine
    engine = FlukaEngine.instance
    if not engine._flukaio_connected:
        raise ValueError(f"Fluka server not yet running!\n"
                        + "Please do this first, by calling xcoll.FlukaEngine.start_server(fluka_input_file.inp).")

    if not engine.has_particle_ref:
        raise ValueError(f"Fluka reference particle not set!\n"
                        + "Please do this first, by calling xcoll.FlukaEngine().set_particle_ref().")

    if 1.4*npart > engine.n_alloc:
        raise ValueError(f"Tracking {npart} particles but only {engine.n_alloc} allocated! "
                       + f"Remember to leaver room for secondaries...")

    engine.init_tracking(npart)

    if particles.particle_id.max() > engine.max_particle_id:
        raise ValueError(f"Some particles have an id that is higher than the highest id known "
                       + f"to FLUKA ({engine.max_particle_id}). This could happen if this "
                       + f"particles object is larger than the the first particles instance "
                       + f"tracked in this session, or if secondary particles are generated "
                       + f"somewhere else than FLUKA. In that case, call "
                       + f"xcoll.FlukaEngine().init_tracking(max_particle_id) before tracking "
                       + f"with a value large enough to accomodate secondaries outside of FLUKA.\n"
                       + f"In any case, please stop and restart the FLUKA server now.")

    particle_ref = engine.particle_ref
    if abs(particles.mass0 - particle_ref.mass0) > 1e-3:
        raise ValueError("Error in reference mass of `particles`: not in sync with reference particle!")
    if abs(particles.q0 - particle_ref.q0) > 1e-3:
        raise ValueError("Error in reference charge of `particles`: not in sync with reference particle!")
    if np.any([pdg_id == 0 for pdg_id in particles.pdg_id]):
        raise ValueError("Some particles are missing the pdg_id!")

    # FLUKA collimators are centered; need to shift
    dx = 0
    dy = 0
    if abs(np.mod(collimator.angle-90,180)-90) < 1e-6:
        dx = collimator.ref_x
    elif abs(np.mod(collimator.angle,180)-90) < 1e-6:
        dy = collimator.ref_y
    # TODO: Let's forget about the closed orbit for the skew collimators for now...
    particles.x -= dx
    particles.y -= dy
    track_core(collimator, particles)
    particles.x += dx
    particles.y += dy


def _expand(arr, dtype=float):
    from .engine import FlukaEngine
    engine = FlukaEngine.instance
    max_part = engine.n_alloc
    return np.concatenate((arr, np.zeros(max_part-arr.size, dtype=dtype)))


def track_core(collimator, part):
    npart = part._num_active_particles
    from .engine import FlukaEngine
    engine = FlukaEngine.instance
    try:
        from .pyflukaf import track_fluka
    except ImportError as error:
        engine._warn_pyfluka(error)
        return

    max_part       = engine.n_alloc
    alive_at_entry = part.state > 0
    max_id         = part.particle_id[alive_at_entry].max()
    assert alive_at_entry.sum() == npart

    # Get particle data
    mass       = part.mass0*part.charge_ratio[alive_at_entry] / part.chi[alive_at_entry]
    charge     = part.q0*part.charge_ratio[alive_at_entry]
    pdg_id     = part.pdg_id[alive_at_entry].copy()
    _, A, Z, _ = pdg.get_properties_from_pdg_id(pdg_id)

    # Prepare arrays for FORTRAN
    data = {}
    data['x']      = _expand(part.x[alive_at_entry].copy() * 1000.)
    data['xp']     = _expand(part.px[alive_at_entry] * part.rpp[alive_at_entry] * 1000.)
    data['y']      = _expand(part.y[alive_at_entry].copy() * 1000.)
    data['yp']     = _expand(part.py[alive_at_entry] * part.rpp[alive_at_entry] * 1000.)
    data['zeta']   = _expand(part.zeta[alive_at_entry].copy() * 1000.)
    data['e']      = _expand(part.energy[alive_at_entry].copy() / 1.e6)
    data['m']      = _expand(mass / 1.e6)
    data['q']      = _expand(charge.astype(np.int16), dtype=np.int16)
    data['A']      = _expand(A.astype(np.int32), dtype=np.int32)
    data['Z']      = _expand(Z.astype(np.int32), dtype=np.int32)
    data['pdg_id'] = _expand(pdg_id.astype(np.int32), dtype=np.int32)
    # FLUKA is 1-indexed
    data['pid']    = _expand(part.particle_id[alive_at_entry].copy().astype(np.int32) + 1, dtype=np.int32)
    data['ppid']   = _expand(part.parent_particle_id[alive_at_entry].copy().astype(np.int32) + 1, dtype=np.int32)
    data['weight'] = _expand(part.weight[alive_at_entry].copy())
    # Hard-coded spin (not used)
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
    track_fluka(turn=turn_in+1,                # Turn indexing start from 1 with FLUKA IO (start from 0 with xpart)
                fluka_id=collimator.fluka_id,
                length=collimator.length,      # FLUKA uses inactive front and back, so pass full length
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

    # IMPORTANT: we assume that a parent can never continue existing after an interaction. TODO: is this correct?

#     # FLUKA does not use a parent ID, but a primary ID (hence not the direct parent but the first impact)
#     # ===================================================================================================
#     # Select particles whose ppid was not sent as input to the collimator
#     # (the parent is from a previous interaction)
#     to_correct = [ppid not in part.particle_id[alive_at_entry] for ppid in new_ppid]
#     assert np.all([ppid in part.parent_particle_id[alive_at_entry] for ppid in new_ppid[to_correct]])

#     # FLUKA returns particles that have undergone a nuclear interaction as new particles: need to correct
#     # ===================================================================================================
#     new_pdg_id = data['pdg_id'][:npart]

#     # We select only particles that have a different pid than their parent, while ensuring the
#     # parent existed as a particle in the input sent to FLUKA (this is to avoid selecting cases
#     # where the input particle already had a different parent from a previous interaction)
#     existed  = [ppid in part.particle_id[alive_at_entry] for ppid in new_ppid]  # TODO: this is slow
#     diff_ids = new_pid != new_ppid

#     # count the number of children per parent
#     _, inv, count_pp = np.unique(new_ppid, return_inverse=True, return_counts=True)
#     only_one_child   = count_pp[inv] == 1

#     # Select only those where the particle type did not change
#     # TODO: sloooooooow
#     old_pdg_id = []
#     for ppid in new_ppid:
#         pdg_id = part.pdg_id[alive_at_entry][part.particle_id[alive_at_entry]==ppid]
#         if len(pdg_id) > 0:
#             old_pdg_id += list(pdg_id[:1])
#         else:
#             old_pdg_id += [0]
#     old_pdg_id = np.array(old_pdg_id)
#     pdg_id_same = new_pdg_id == old_pdg_id

#     # re-assign the original particle ids if the conditions are met
#     nucl_int = existed & diff_ids & only_one_child & pdg_id_same
#     new_pid[nucl_int] = new_ppid[nucl_int]


    # TODO: Impact Table
    # Absorbed: trivial
    # Not hit:    pid+ppid did not change and dE = 0
    # Hit (MCS):  pid+ppid did not change and dE != 0
    # Hit (nucl): above mask + children


    # Update existing particles  (these missed the collimator or only underwent elastic interactions)
    # ===============================================================================================
    mask_new = new_pid <= max_id

    if np.any(mask_new):
        idx_old  = np.array([np.where(part.particle_id==idx)[0][0] for idx in new_pid[mask_new]])  # list of indices

        # Sanity check
        assert np.all(part.particle_id[idx_old] == new_pid[mask_new])
        assert np.all(part.parent_particle_id[idx_old] == new_ppid[mask_new])
        assert np.all(part.state[idx_old] > 0)

        # Update energy   TODO: in principle FLUKA uses pc, not energy!!
        E_diff = np.zeros(len(part.x))
        E_diff[idx_old] = data['e'][:npart][mask_new]*1.e6 - part.energy[idx_old]
        part.add_to_energy(E_diff) # TODO: need to correct for weight
        collimator.accumulated_energy -= E_diff.sum()
        rpp = part.rpp[idx_old]    # This is now already updated by the new energy

        # Update other fields
        part.x[idx_old]            = data['x'][:npart][mask_new] / 1000.
        part.px[idx_old]           = data['xp'][:npart][mask_new] / rpp / 1000.
        part.y[idx_old]            = data['y'][:npart][mask_new] / 1000.
        part.py[idx_old]           = data['yp'][:npart][mask_new] / rpp / 1000.
        part.zeta[idx_old]         = data['zeta'][:npart][mask_new] / 1000.
        part.charge_ratio[idx_old] = data['q'][:npart][mask_new] / part.q0
        part.chi[idx_old]          = data['q'][:npart][mask_new] / (data['m'][:npart][mask_new]*1.e6) \
                                               * part.mass0 / part.q0
        part.pdg_id[idx_old]       = data['pdg_id'][:npart][mask_new]
        part.weight[idx_old]       = data['weight'][:npart][mask_new]
        part.s[idx_old]            = s_in + collimator.length

    # Little hack to set the dead particles, as idx_old is not a mask (but a list of indices)
    # (hence we cannot use ~idx_old)
    part.state[alive_at_entry] = -335  # XC_LOST_ON_FLUKA
    if np.any(mask_new):
        part.state[idx_old]        = 1     # These actually survived
    collimator.accumulated_energy += part.energy[alive_at_entry & (part.state==-335)].sum() # TODO: need to correct for weight

#     # Correct state for parent particles that disappeared because of multiple children
#     mask = [pid not in new_pid and pid in new_ppid for pid in part.particle_id]
#     part.state[alive_at_entry & mask] = -385  # XC_PARENT_ON_FLUKA

#     # Sanity check
#     was_alive = [ppid in part.particle_id[alive_at_entry] for ppid in new_ppid]
#     _, inv, count_pp = np.unique(new_ppid, return_inverse=True, return_counts=True)
#     more_than_one_child = count_pp[inv] > 1
#     assert set(part.particle_id[alive_at_entry & mask]) == set(new_ppid[was_alive & more_than_one_child])


    # Add new particles
    # ================
    mask_new = new_pid > max_id

    if not np.any(mask_new):
#         # Sanity check
#         assert -385 not in part.state[alive_at_entry]
        part.reorganize()

    else:
#         # Sanity check
#         print(new_ppid[mask_new])
#         print(new_pid[mask_new])
#         for pid, ppid in zip(new_pid[mask_new], new_ppid[mask_new]):
#             if ppid not in part.particle_id[alive_at_entry]:
#                 mask = part.particle_id == ppid
#                 print(f"Oeps: {ppid}->{pid} not found. State {part.state[mask]}.")
#         assert np.all([ppid in part.particle_id[alive_at_entry] for ppid in new_ppid[mask_new]])
#         assert np.all([part.state[part.particle_id==ppid][0] == -385 for ppid in new_ppid[mask_new]])

        # Check that there is enough room in the particles object
        num_free = part._capacity - part._num_lost_particles - part._num_active_particles
        num_needed = np.sum(mask_new)
        if num_free < num_needed:
            extra_capacity = int(1.2*(num_needed - num_free))  # 20% margin
            # TODO: this does not work!!
            part = xp.Particles.from_dict(part.to_dict(), _capacity=part._capacity+extra_capacity)

        # Create new particles
        # TODO: FLUKA uses pc; then calculate delta from p
        m     = data['m'][:npart][mask_new] * 1.e6
        E     = data['e'][:npart][mask_new] * 1.e6
        E0    = part.energy0[0]
        b0    = part.beta0[0]
        m0    = part.mass0
        delta = np.sqrt(1./(b0**2) * E**2/(E0**2) * m0**2/(m**2) - 1./(b0**2) + 1.) - 1.
        rpp   = 1. / (1. + delta)
        new_particles = xp.Particles(_context=part._buffer.context,
                p0c = part.p0c[0],
                mass0 = part.mass0,
                q0 = part.q0,
                s = s_in + collimator.length,
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
#                 particle_id = new_pid[mask_new],  # do not use the FLUKA pid as it will keep on increasing more
                parent_particle_id = new_ppid[mask_new],
                weight = data['weight'][:npart][mask_new],
                start_tracking_at_element = start)

# TODO: energy of parents is wrong !
#         # Correct energy transfer: parent does not have energy of children anymore
#         idx_old = np.array([np.where(part.particle_id==idx)[0][0] for idx in new_particles.parent_particle_id])
#         E_diff = np.zeros(len(part.x))
#         np.add.at(E_diff, idx_old, -new_particles.energy)
#         part.add_to_energy(E_diff)
        # This energy was added too much to the accumulation (by the -335 state), so subtract (E_diff < 0)
        collimator.accumulated_energy -= new_particles.energy.sum() # TODO: need to correct for weight

        # Add new particles
        new_particles._init_random_number_generator()
        part.add_particles(new_particles)
        engine.max_particle_id = part.particle_id.max()







