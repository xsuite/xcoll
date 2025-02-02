# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pytest

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc
import time

import matplotlib.pyplot as plt
from scipy.stats import ks_2samp


@pytest.mark.parametrize('num_part', [1000, 5000])
def test_simple_track(num_part):
    print(f"Running test_simple_track with {num_part} particles")
    _capacity = num_part*2

    # If a previous test failed, stop the server manually
    if xc.FlukaEngine.is_running():
        xc.FlukaEngine.stop(clean=True)

    # Define collimator and start the FLUKA server
    coll = xc.FlukaCollimator(length=0.6, jaw=0.001, assembly='hilumi_tcppm')
    coll_name = 'tcp.c6l7.b1'
    xc.FlukaEngine.start(elements=coll, names=coll_name, debug_level=1, _capacity=_capacity)

    # Particle distribution
    x_init   = np.random.normal(loc=0.002, scale=1e-3, size=num_part)
    px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    xc.FlukaEngine.set_particle_ref(particle_ref=particle_ref)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                   particle_ref=particle_ref, _capacity=_capacity)
    part_fluka = part_init.copy()
    part_drift = part_init.copy()

    # FLUKA tracking
    start = time.time()
    coll.track(part_fluka)
    fluka_time = round(time.time()-start, 3)
    expected_time = round(0.025 * num_part)
    if fluka_time > expected_time:
        raise RuntimeError(f"FLUKA tracking took {fluka_time}s, exceeding threshold of {expected_time}s")
    else:
        print(f"FLUKA tracking took {fluka_time}s")

    # Drift tracking
    coll._equivalent_drift.length = coll.length
    coll._equivalent_drift.track(part_drift)

    mask_fluka = part_fluka.state > 0
    mask_drift = part_drift.state > 0
    ks_stat, p_value = ks_2samp(part_fluka.x[mask_fluka], part_drift.x[mask_drift])
    assert p_value <= 0.05, f"Distributions are not significantly different (p = {p_value})"
    print(f"KS test passed with p = {p_value}")
    mask_fluka = part_fluka.state > 0
    mask_drift = part_drift.state > 0
    ks_stat, p_value = ks_2samp(part_fluka.x[mask_fluka], part_drift.x[mask_drift])
    assert p_value <= 0.05, f"Distributions are not significantly different (p = {p_value})"
    print(f"KS test passed with p = {p_value}")

    # Stop the FLUKA server
    xc.FlukaEngine.stop(clean=True)


# @pytest.mark.parametrize('tilt', [0, [2.2e-6, 1.3e-6], [1.9e-6, -2.7e-6]],
#                          ids=['no_tilt', 'positive_tilt', 'pos_neg_tilt'])
@pytest.mark.parametrize('assembly', ['sps_tcsm', 'lhc_tcp', 'lhc_tcpm', 'lhc_tcsg',
                                      'lhc_tcsp', 'lhc_tcla', 'lhc_tct', 'lhc_tcl',
                                      'lhc_tclia', 'lhc_tclib', 'lhc_tcdqaa', 'lhc_tcdqab',
                                      'lhc_tcdqac', 'hilumi_tcld', 'hilumi_tctx', 'hilumi_tcty',
                                      'hilumi_tctpx', 'hilumi_tclx', 'fcc_tcp', 'fcc_tcsg'])
@pytest.mark.parametrize('angle', [0, 90, 130.5])
@pytest.mark.parametrize('jaw', [0.001, [0.0013, -0.002789], [-1.2e-6, -3.2e-3], [3.789e-3, 4.678e-7]],
                         ids=['symmetric', 'asymmetric', 'negative', 'positive'])
def test_fluka_jaw(jaw, angle, assembly):
    _JAW_ACCURACY = 1e-12  # Anything in this region around the jaw might or might not hit; we can't be sure
    _MOMENTUM_ACCURACY = 1e-12
    num_part = 5000
    _capacity = num_part*2
    jaw_band = 1.e-6
    tilt = 0

    # If a previous test failed, stop the server manually
    if xc.FlukaEngine.is_running():
        xc.FlukaEngine.stop(clean=True)

    # Define collimator and start the FLUKA server
    side = 'left' if 'tcdq' in assembly else 'both'
    coll = xc.FlukaCollimator(length=1, jaw=jaw, angle=angle, tilt=tilt, assembly=assembly, side=side)
    coll_name = 'tcp.c6l7.b1'
    xc.FlukaEngine.start(elements=coll, names=coll_name, debug_level=1, _capacity=_capacity)
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    xc.FlukaEngine.set_particle_ref(particle_ref=particle_ref)

    # Particle distribution (x and y are in the frame of the collimator)
    num_part_step = num_part//5
    num_part = 5*num_part_step
    x = np.random.uniform(-0.02, 0.02, num_part_step)
    if coll.side != 'both':
        num_part_step *= 2
    if coll.side == 'left' or coll.side == 'both':
        jaw_L = min(coll.jaw_LU, coll.jaw_LD)
        x = np.concatenate([x, np.random.uniform(jaw_L - jaw_band, jaw_L -_JAW_ACCURACY, num_part_step)])
        x = np.concatenate([x, np.random.uniform(jaw_L +_JAW_ACCURACY, jaw_L + jaw_band, num_part_step)])
    else:
        jaw_L = 1e6
    if coll.side == 'right' or coll.side == 'both':
        jaw_R = max(coll.jaw_RU, coll.jaw_RD)
        x = np.concatenate([x, np.random.uniform(jaw_R - jaw_band, jaw_R -_JAW_ACCURACY, num_part_step)])
        x = np.concatenate([x, np.random.uniform(jaw_R +_JAW_ACCURACY, jaw_R + jaw_band, num_part_step)])
    else:
        jaw_R = -1e6
    y = np.random.uniform(-0.02, 0.02, num_part)
    x_new = np.cos(np.deg2rad(angle))*x - np.sin(np.deg2rad(angle))*y
    y_new = np.sin(np.deg2rad(angle))*x + np.cos(np.deg2rad(angle))*y
    part_init = xp.build_particles(x=x_new, y=y_new, particle_ref=xc.FlukaEngine().particle_ref,
                                   _capacity=xc.FlukaEngine()._capacity)
    mask = np.concatenate([(x >= jaw_L) | (x <= jaw_R),
                          np.full(_capacity-num_part, False)])
    hit_ids = part_init.particle_id[mask & (part_init.state > 0)]
    not_hit_ids = part_init.particle_id[~mask & (part_init.state > 0)]

    # TODO: jaw tilts, and particle angles

    # Track
    part = part_init.copy()
    coll.track(part)
    mask_hit = np.isin(part.parent_particle_id, hit_ids)
    mask_not_hit = np.isin(part.parent_particle_id, not_hit_ids)

    # mask = mask_not_hit & (part.state < 1)
    # wrong_ids = part.parent_particle_id[mask]
    # mask_wrong_init = np.isin(part_init.particle_id, wrong_ids)

    # mask_hit_init = np.isin(part_init.particle_id, hit_ids)
    # mask_not_hit_init = np.isin(part_init.particle_id, not_hit_ids)
    # plt.scatter(part_init.x[mask_not_hit_init], part_init.y[mask_not_hit_init], c='b', s=2)
    # plt.scatter(part_init.x[mask_hit_init], part_init.y[mask_hit_init], c='g', s=2)
    # plt.scatter(part_init.x[mask_wrong_init], part_init.y[mask_wrong_init], c='r', s=2)
    # plt.axvline(coll.jaw_LU, c='k', linestyle='--')
    # plt.axvline(coll.jaw_LD, c='k', linestyle='--')
    # plt.axvline(coll.jaw_RU, c='k', linestyle='--')
    # plt.axvline(coll.jaw_RD, c='k', linestyle='--')
    # plt.show()


    # Particles that are supposed to not have hit the collimator, but have a kick or are dead, are considered faulty
    assert not np.any(abs(part.px[mask_not_hit]) > _MOMENTUM_ACCURACY)
    assert not np.any(abs(part.py[mask_not_hit]) > _MOMENTUM_ACCURACY)
    assert not np.any(part.state[mask_not_hit] < 1)

    # Particles that are supposed to have hit the collimator, but are alive and have no kick, are considered faulty
    faulty =  mask_hit & (abs(part.px) < _MOMENTUM_ACCURACY) & (abs(part.py) < _MOMENTUM_ACCURACY)
    faulty &= (part.state > 0)
    assert len(part.x[faulty]) <= 1  # We allow for a small margin of error

    # Stop the FLUKA server
    xc.FlukaEngine.stop(clean=True)
