# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import pytest

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc
import time

import matplotlib.pyplot as plt

# TODO: particle angles

jaws = [0.001, [0.0013, -0.002789], [-1.2e-6, -3.2e-3], [3.789e-3, 4.678e-7]]
jaw_ids = ['symmetric', 'asymmetric', 'negative', 'positive']
angles = [0, 90, 127.5]
tilts = [0, [2.2e-6, 1.3e-6], [1.9e-6, -2.7e-6]]
tilt_ids = ['no_tilt', 'positive_tilt', 'pos_neg_tilt']

particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)


@pytest.mark.parametrize('tilt', tilts, ids=tilt_ids)
@pytest.mark.parametrize('angle', angles)
@pytest.mark.parametrize('jaw', jaws, ids=jaw_ids)
def test_everest(jaw, angle, tilt):
    coll = xc.EverestCollimator(length=1, jaw=jaw, angle=angle, tilt=tilt, material=xc.materials.MolybdenumGraphite)
    part_init, hit_ids, not_hit_ids = _generate_particles(coll, num_part=250_000, particle_ref=particle_ref)
    part = part_init.copy()
    coll.track(part)
    _assert_valid_positions(part, hit_ids, not_hit_ids)


# TODO: lhc_tdi and fcc_tcdq still fail. Are they unrotatable? Side fixed or ill-defined?
# @pytest.mark.parametrize('tilt', tilts, ids=tilt_ids)
@pytest.mark.parametrize('assembly', ['sps_tcsm', 'lhc_tcp', 'lhc_tcsg', 'lhc_tcsp',
                                      'lhc_tcla', 'lhc_tct', 'lhc_tcl', 'lhc_tdi',
                                      'lhc_tclia', 'lhc_tclib', 'lhc_tcdqaa', 'lhc_tcdqab',
                                      'lhc_tcdqac', 'hilumi_tcppm', 'hilumi_tcspm', 'hilumi_tcsg',
                                      'hilumi_tcld', 'hilumi_tctx', 'hilumi_tcty', 'hilumi_tclx',
                                      'fcc_tcp', 'fcc_tcsg', 'fcc_tcdq'])
@pytest.mark.parametrize('angle', angles)
@pytest.mark.parametrize('jaw', jaws, ids=jaw_ids)
def test_fluka(jaw, angle, assembly):
    tilt = 0

    # If a previous test failed, stop the server manually
    if xc.FlukaEngine.is_running():
        xc.FlukaEngine.stop(clean=True)

    # Define collimator and start the FLUKA server
    coll = xc.FlukaCollimator(length=1, jaw=jaw, angle=angle, tilt=tilt, assembly=assembly)
    xc.FlukaEngine.particle_ref = particle_ref
    xc.FlukaEngine.start(elements=coll, capacity=10_000)

    part_init, hit_ids, not_hit_ids = _generate_particles(coll, num_part=5000, dim=0.015,
                            _capacity=xc.FlukaEngine.capacity, particle_ref=xc.FlukaEngine().particle_ref)
    part = part_init.copy()
    coll.track(part)
    _assert_valid_positions(part, hit_ids, not_hit_ids)

    # Stop the FLUKA server
    xc.FlukaEngine.stop(clean=True)


def _generate_particles(coll, num_part, particle_ref, _capacity=None, jaw_band=1.e-6,
                        jaw_accuracy=1.e-12, dim=0.05):
    if _capacity is None:
        _capacity = num_part

    # Particle distribution (x and y are in the frame of the collimator)
    num_part_step = num_part//5
    num_part = 5*num_part_step
    x = np.random.uniform(-dim, dim, num_part_step)
    if coll.side != 'both':
        num_part_step *= 2
    if coll.side == 'left' or coll.side == 'both':
        jaw_L = min(coll.jaw_LU, coll.jaw_LD)
        x = np.concatenate([x, np.random.uniform(jaw_L - jaw_band, jaw_L - jaw_accuracy, num_part_step)])
        x = np.concatenate([x, np.random.uniform(jaw_L + jaw_accuracy, jaw_L + jaw_band, num_part_step)])
    else:
        jaw_L = 1e6
    if coll.side == 'right' or coll.side == 'both':
        jaw_R = max(coll.jaw_RU, coll.jaw_RD)
        x = np.concatenate([x, np.random.uniform(jaw_R - jaw_band, jaw_R - jaw_accuracy, num_part_step)])
        x = np.concatenate([x, np.random.uniform(jaw_R + jaw_accuracy, jaw_R + jaw_band, num_part_step)])
    else:
        jaw_R = -1e6
    y = np.random.uniform(-dim, dim, num_part)
    x_new = np.cos(np.deg2rad(coll.angle))*x - np.sin(np.deg2rad(coll.angle))*y
    y_new = np.sin(np.deg2rad(coll.angle))*x + np.cos(np.deg2rad(coll.angle))*y
    part_init = xp.build_particles(x=x_new, y=y_new, particle_ref=particle_ref, _capacity=_capacity)
    mask = np.concatenate([(x >= jaw_L) | (x <= jaw_R),
                          np.full(_capacity-num_part, False)])
    hit_ids = part_init.particle_id[mask & (part_init.state > 0)]
    not_hit_ids = part_init.particle_id[~mask & (part_init.state > 0)]

    return part_init, hit_ids, not_hit_ids


def _assert_valid_positions(part, hit_ids, not_hit_ids, momentum_accuracy=1.e-12):
    mask_hit = np.isin(part.parent_particle_id, hit_ids)
    mask_not_hit = np.isin(part.parent_particle_id, not_hit_ids)

    # Particles that are supposed to not have hit the collimator, but have a kick or are dead, are considered faulty
    assert not np.any(abs(part.px[mask_not_hit]) > momentum_accuracy)
    assert not np.any(abs(part.py[mask_not_hit]) > momentum_accuracy)
    assert not np.any(part.state[mask_not_hit] < 1)

    # Particles that are supposed to have hit the collimator, but are alive and have no kick, are considered faulty
    faulty =  mask_hit & (abs(part.px) < momentum_accuracy) & (abs(part.py) < momentum_accuracy)
    faulty &= (part.state > 0)
    assert len(part.x[faulty]) <= 1  # We allow for a small margin of error


def _plot_jaws(coll, part_init, part, hit_ids, not_hit_ids):
    mask = np.isin(part.parent_particle_id, not_hit_ids) & (part.state < 1)
    wrong_ids = part.parent_particle_id[mask]
    mask_wrong_init = np.isin(part_init.particle_id, wrong_ids)

    mask_hit_init = np.isin(part_init.particle_id, hit_ids)
    mask_not_hit_init = np.isin(part_init.particle_id, not_hit_ids)
    plt.scatter(part_init.x[mask_not_hit_init], part_init.y[mask_not_hit_init], c='b', s=2)
    plt.scatter(part_init.x[mask_hit_init], part_init.y[mask_hit_init], c='g', s=2)
    plt.scatter(part_init.x[mask_wrong_init], part_init.y[mask_wrong_init], c='r', s=2)
    plt.axvline(coll.jaw_LU, c='k', linestyle='--')
    plt.axvline(coll.jaw_LD, c='k', linestyle='--')
    plt.axvline(coll.jaw_RU, c='k', linestyle='--')
    plt.axvline(coll.jaw_RD, c='k', linestyle='--')
    plt.show()
