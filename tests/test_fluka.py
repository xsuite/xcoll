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
    xc.FlukaEngine.particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    xc.FlukaEngine.start(elements=coll, names=coll_name, capacity=_capacity, verbose=True)

    # Particle distribution
    x_init   = np.random.normal(loc=0.002, scale=1e-3, size=num_part)
    px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                   particle_ref=xc.FlukaEngine.particle_ref,
                                   _capacity=xc.FlukaEngine.capacity)
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


def test_prototypes():
    raise ValueError("Need to write test for FlukaAssembly to check registry works as expected")

