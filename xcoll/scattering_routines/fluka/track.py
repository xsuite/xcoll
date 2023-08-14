# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

from ...impacts import make_ion_from_properties

def drift_4d(x, y, xp, yp, length):
    x += xp * length
    y += yp * length
    return

def drift_zeta(zeta, rvv, xp, yp, length):
    rv0v = 1./rvv
    dzeta = 1 - rv0v * (1 + (xp**2 + yp**2)/2 )
    zeta += length * dzeta
    return zeta

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

        from .engine import FlukaEngine
        engine = FlukaEngine.instance
        if not engine._flukaio_connected:
            raise ValueError(f"Fluka server not yet running!\n"
                            + "Please do this first, by calling xcoll.FlukaEngine.start_server(fluka_input_file.inp).")
        npart = particles._num_active_particles
        if npart > FlukaEngine.instance.n_alloc:
            raise ValueError(f"Tracking {npart} particles but only {FlukaEngine.instance.n_alloc} allocated!")
        if npart == 0:
            return

        # TODO: when in C, drifting should call routine from drift element
        #       such that can be exact etc
        if not collimator.active or not collimator._tracking:
            # Drift full length
            drift_6d(particles, collimator.length)
        else:
            # FLUKA uses inactive front and back, so pass full length
            track_core(collimator, collimator.length)
            particles.reorganize()
        return


def track_core(collimator, part):
    npart = part._num_active_particles
    from .engine import FlukaEngine
    engine = FlukaEngine.instance
    try:
        from .pyflukaf import track_fluka
    except ImportError as error:
        engine._warn_pyfluka(error)
        return

    # Add a margin of three times the particles size for secondaries
    max_part = 6*npart

    # Get particle data
    data = {}
    data['x_part']      = np.concatenate((part.x[:npart].copy() * 1000.,              np.zeros(max_part-npart, dtype=float)))
    data['xp_part']     = np.concatenate((part.px[:npart] * part.rpp[:npart] * 1000., np.zeros(max_part-npart, dtype=float)))
    data['y_part']      = np.concatenate((part.y[:npart].copy() * 1000.,              np.zeros(max_part-npart, dtype=float)))
    data['yp_part']     = np.concatenate((part.py[:npart] * part.rpp[:npart] * 1000., np.zeros(max_part-npart, dtype=float)))
    data['zeta_part']   = np.concatenate((part.zeta[:npart].copy() * 1000.,           np.zeros(max_part-npart, dtype=float)))
    data['e_part']      = np.concatenate((part.energy[:npart].copy() / 1.e6,          np.zeros(max_part-npart, dtype=float)))
    mass   = part.mass0*part.charge_ratio[:npart]/part.chi[:npart] / 1.e6
    charge = part.q0*part.charge_ratio[:npart]
    A, Z, pdgid = make_ion_from_properties(charge, mass)
    data['m_part']      = np.concatenate((mass,                                      np.zeros(max_part-npart, dtype=float)))
    data['q_part']      = np.concatenate((charge.astype(np.int32),                   np.zeros(max_part-npart, dtype=np.int32)))
    data['A_part']      = np.concatenate((A,                                         np.zeros(max_part-npart, dtype=np.int32)))
    data['Z_part']      = np.concatenate((Z,                                         np.zeros(max_part-npart, dtype=np.int32)))
    data['pdgid_part']  = np.concatenate((pdgid,                                     np.zeros(max_part-npart, dtype=np.int32)))
    # TODO: hard-coded spin
    data['spin_x_part'] = np.concatenate((np.ones(npart)*0.5,                        np.zeros(max_part-npart, dtype=float)))
    data['spin_y_part'] = np.concatenate((np.ones(npart)*0.5,                        np.zeros(max_part-npart, dtype=float)))
    data['spin_z_part'] = np.concatenate((np.ones(npart)*0.5,                        np.zeros(max_part-npart, dtype=float)))
    data['partID']      = np.concatenate((part.particle_id[:npart].copy().astype(np.int32),
                                                                                     np.zeros(max_part-npart, dtype=np.int32)))
    data['parentID']    = np.concatenate((part.parent_particle_id[:npart].copy().astype(np.int32),
                                                                                     np.zeros(max_part-npart, dtype=np.int32)))
    data['partWeight']  = np.concatenate((part.weight[:npart].copy(),                np.zeros(max_part-npart, dtype=float)))
    
    # send to fluka
    track_fluka(turn=part.at_turn[0],
                fluka_id=collimator.collimator_id,
                length=collimator.length,
                alive_part=npart, # max_part,
                x_part=data['x_part'],
                xp_part=data['xp_part'],
                y_part=data['y_part'],
                yp_part=data['yp_part'],
                zeta_part=data['zeta_part'],
                e_part=data['e_part'],
                m_part=data['m_part'],
                q_part=data['q_part'],
                a_part=data['A_part'],
                z_part=data['Z_part'],
                pdgid_part=data['pdgid_part'],
                part_id=data['partID'],
                parent_id=data['parentID'],
                part_weight=data['partWeight'],
                spin_x_part=data['spin_x_part'],
                spin_y_part=data['spin_y_part'],
                spin_z_part=data['spin_z_part']
               )

    return  npart, max_part, \
            data['x_part'], data['xp_part'], data['y_part'], data['yp_part'], data['zeta_part'],\
            data['e_part'], data['m_part'], data['q_part'], data['A_part'], data['Z_part'],\
            data['pdgid_part'], data['partID'], data['parentID'], data['partWeight'],\
            data['spin_x_part'], data['spin_y_part'], data['spin_z_part']







