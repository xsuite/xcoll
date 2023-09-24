# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np
import xpart as xp
import xpart.pdg as pdg


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

    max_part = engine.n_alloc
    max_id   = part.particle_id.max()

    # Get particle data
    mass       = part.mass0*part.charge_ratio[:npart]/part.chi[:npart]
    charge     = part.q0*part.charge_ratio[:npart]
    pdg_id     = part.pdg_id[:npart].copy()
    _, A, Z, _ = pdg.get_properties_from_pdg_id(pdg_id)

    # Prepare arrays for FORTRAN
    data = {}
    data['x']      = _expand(part.x[:npart].copy() * 1000.)
    data['xp']     = _expand(part.px[:npart] * part.rpp[:npart] * 1000.)
    data['y']      = _expand(part.y[:npart].copy() * 1000.)
    data['yp']     = _expand(part.py[:npart] * part.rpp[:npart] * 1000.)
    data['zeta']   = _expand(part.zeta[:npart].copy() * 1000.)
    data['e']      = _expand(part.energy[:npart].copy() / 1.e6)
    data['m']      = _expand(mass / 1.e6)
    data['q']      = _expand(charge.astype(np.int32), dtype=np.int32)
    data['A']      = _expand(A.astype(np.int32), dtype=np.int32)
    data['Z']      = _expand(Z.astype(np.int32), dtype=np.int32)
    data['pdg_id'] = _expand(pdg_id.astype(np.int32), dtype=np.int32)
    data['pid']    = _expand(part.particle_id[:npart].copy().astype(np.int32)+1, dtype=np.int32)        # FLUKA is 1-indexed
    data['ppid']   = _expand(part.parent_particle_id[:npart].copy().astype(np.int32)+1, dtype=np.int32) # FLUKA is 1-indexed
    data['weight'] = _expand(part.weight[:npart].copy())
    # Hard-coded spin (not used)
    data['spin_x'] = _expand(np.zeros(npart))
    data['spin_y'] = _expand(np.zeros(npart))
    data['spin_z'] = _expand(np.zeros(npart))

    # Change npart to np.array to make it writable, store some initial data
    npart_in = npart
    npart    = np.array(npart, dtype=np.int64)
    s_in     = part.s[0]
    ele_in   = part.at_element[0]
    turn_in  = part.at_turn[0]

    # send to fluka
    track_fluka(turn=part.at_turn[0]+1,       # Turn indexing start from 1 with FLUKA IO (start from 0 with xpart)
                fluka_id=collimator.fluka_id,
                length=collimator.length,     # FLUKA uses inactive front and back, so pass full length
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
                pdgid_part=data['pdg_id'],
                part_id=data['pid'],
                parent_id=data['ppid'],
                part_weight=data['weight'],
                spin_x_part=data['spin_x'],
                spin_y_part=data['spin_y'],
                spin_z_part=data['spin_z']
               )

    new_pid    = data['pid'][:npart] - 1   # return to python 0-index
    new_ppid   = data['ppid'][:npart] - 1  # return to python 0-index
    new_pdg_id = data['pdg_id'][:npart]

    idx_old    = np.array([np.where(part.particle_id[:npart_in]==idx)[0][0] for idx in new_ppid]) # TODO: clean numpy-style?
    old_pdg_id = part.pdg_id[:npart_in][idx_old]

    # FLUKA returns particles that have undergone a nuclear interaction as new particles
    # We don't want that. Hence, when a particle id is different from its parent id, but
    # the pdg_id did not change and no other children are created, we reset the particle
    # id to the parent id
    # TODO: is this robust????? What if the child is absorbed in the same collimator as it was created in????
    diff_ids = new_pid != new_ppid  # particles with different ids
    # count the number of children per parent
    _, inv, count_pp = np.unique(new_ppid, return_inverse=True, return_counts=True)
    only_one_child   = count_pp[inv] == 1
    # Select only those where the particle type did not change
    pdg_id_same = new_pdg_id == old_pdg_id
    # re-assign the original particle ids if the conditions are met
    mask_pid = diff_ids & only_one_child & pdg_id_same
    new_pid[mask_pid] = new_ppid[mask_pid]

    # TODO:
    # How to know when hit?
    # When energy before/after is different

    # Update existing particles
    # =========================
    mask_new = new_pid <= max_id
    idx_old  = np.array([np.where(part.particle_id==idx)[0][0] for idx in new_pid[mask_new]])

    # Sanity check
    assert np.all(part.particle_id[:npart_in][idx_old] == new_pid[mask_new])
    assert np.all(part.parent_particle_id[:npart_in][idx_old] == new_ppid[mask_new])

    # Update energy   TODO: this should be nicer, though we needed the automatic update of delta etc
    E_new = dict(zip(new_pid, data['e']*1.e6))
    E_diff = np.array([
            E_new[pid] - E if pid in new_pid[mask_new] else 0.
            for pid, E in zip(part.particle_id, part.energy)
    ])
    part.add_to_energy(E_diff)
    rpp  = part.rpp[:npart_in][idx_old]

    # Update other fields
    part.x[:npart_in][idx_old]            = data['x'][:npart][mask_new] / 1000.
    part.px[:npart_in][idx_old]           = data['xp'][:npart][mask_new] / rpp / 1000.
    part.y[:npart_in][idx_old]            = data['y'][:npart][mask_new] / 1000.
    part.py[:npart_in][idx_old]           = data['yp'][:npart][mask_new] / rpp / 1000.
    part.zeta[:npart_in][idx_old]         = data['zeta'][:npart][mask_new] / 1000.
    part.charge_ratio[:npart_in][idx_old] = data['q'][:npart][mask_new] / part.q0
    part.chi[:npart_in][idx_old]          = data['q'][:npart][mask_new] / (data['m'][:npart][mask_new]*1.e6) \
                                           * part.mass0 / part.q0 
    part.pdg_id[:npart_in][idx_old]       = data['pdg_id'][:npart][mask_new]
    part.weight[:npart_in][idx_old]       = data['weight'][:npart][mask_new]
    part.s[:npart_in][idx_old]            = s_in + collimator.length

    # Little hack to set the dead particles, as idx_old is not a mask (but a list of indices)
    # (hence we cannot use ~idx_old)
    part.state[:npart_in]                 = -335  # XC_LOST_ON_FLUKA
    part.state[:npart_in][idx_old]        = 1     # These actually survived

    # Correct state for parent particles that disappeared because of multiple children
    mask = [pid not in new_pid and pid in new_ppid for pid in part.particle_id[:npart_in]]
    part.state[:npart_in][mask]            = -380  # XC_PARENT_TO_MULTIPLE
    # Sanity check
    _, inv, count_pp = np.unique(new_ppid, return_inverse=True, return_counts=True)
    more_than_one_child = count_pp[inv] > 1
    assert set(part.particle_id[:npart_in][mask]) == set(new_ppid[more_than_one_child])


    # Add new particles
    # ================
    mask_new = new_pid > max_id
    if not np.any(mask_new):
        part.reorganize()
    else:
        rpp = part.rpp[:npart][mask_new]
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
                energy = data['e'][:npart][mask_new] * 1.e6,
                mass_ratio = data['m'][:npart][mask_new]*1.e6 / part.mass0,
                charge_ratio = data['q'][:npart][mask_new] / part.q0,
                at_element = ele_in,
                at_turn = turn_in,
#                 particle_id = new_pid[mask_new],  # do not use the FLUKA pid as it will keep on increasing more
                parent_particle_id = new_ppid[mask_new],
                weight = data['weight'][:npart][mask_new])
        part.add_particles(new_particles)







