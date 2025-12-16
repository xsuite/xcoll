# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import pytest
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc

from xcoll.scattering_routines.fluka.environment import format_fluka_float


@pytest.mark.parametrize('num_part', [1000, 5000])
def test_simple_track(num_part):
    print(f"Running test_simple_track with {num_part} particles")
    _capacity = num_part*2

    # If a previous test failed, stop the server manually
    if xc.fluka.engine.is_running():
        xc.fluka.engine.stop(clean=True)

    # Define collimator and start the FLUKA server
    coll = xc.FlukaCollimator(length=0.6, jaw=0.001, assembly='hilumi_tcppm')
    xc.fluka.engine.particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    xc.fluka.engine.start(elements=coll, capacity=_capacity, clean=False, verbose=True)

    # Particle distribution
    x_init   = np.random.normal(loc=0.002, scale=1e-3, size=num_part)
    px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                   particle_ref=xc.fluka.engine.particle_ref,
                                   _capacity=xc.fluka.engine.capacity)
    part_fluka = part_init.copy()
    part_drift = part_init.copy()

    # FLUKA tracking
    start = time.time()
    coll.track(part_fluka)
    fluka_time = round(time.time()-start, 3)
    # expected_time = round(0.025 * num_part)
    # if fluka_time > expected_time:
    #     raise RuntimeError(f"FLUKA tracking took {fluka_time}s, exceeding threshold of {expected_time}s")
    # else:
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
    xc.fluka.engine.stop(clean=True)


def test_fluka_format_float():
    for i in [0.0000000000011229964381065508,
              0.000000000011229964381065508,
              0.00000000011229964381065508,
              0.0000000011229964381065508,
              0.000000011229964381065508,
              0.00000011229964381065508,
              0.0000011229964381065508,
              0.000011229964381065508,
              0.00011229964381065508,
              0.0011229964381065508,
              0.011229964381065508,
              0.11229964381065508,
              1.1229964381065508,
              11.229964381065508,
              112.29964381065508,
              1122.9964381065508,
              11229.964381065508,
              112299.64381065508,
              1122996.4381065508,
              11229964.381065508,
              112299643.81065508,
              1122996438.1065508,
              11229964381.065508,
              112299643810.65508,
              1122996438106.5508,
              1.12299643810655070652e-12,
              1.12299643810655078730e-11,
              1.12299643810655078730e-10,
              1.12299643810655070975e-09,
              1.12299643810655087519e-08,
              1.12299643810655074284e-07,
              1.12299643810655090166e-06,
              1.12299643810655073225e-05,
              1.12299643810655073225e-04,
              1.12299643810655084067e-03,
              1.12299643810655075393e-02,
              1.12299643810655075393e-01,
              1.12299643810655069842e+00,
              1.12299643810655087606e+01,
              1.12299643810655084053e+02,
              1.12299643810655084053e+03,
              1.12299643810655088600e+04,
              1.12299643810655077687e+05,
              1.12299643810655083507e+06,
              1.12299643810655083507e+07,
              1.12299643810655087233e+08,
              1.12299643810655069351e+09,
              1.12299643810655078888e+10,
              1.12299643810655075073e+11,
              1.12299643810655078125e+12]:
        assert len(format_fluka_float(i)) == 10
        assert len(format_fluka_float(-i)) == 10


def test_particle_ids():
    if xc.fluka.engine.is_running():
        xc.fluka.engine.stop(clean=True)
    coll = xc.FlukaCollimator(length=0.0001, assembly='donadonj', jaw=0)
    xc.fluka.engine.particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='positron', p0c=6.8e12)
    xc.fluka.engine.capacity = 100_000
    xc.fluka.engine.seed = 7856231
    xc.fluka.engine.start(elements=coll, clean=False, verbose=True)
    x_init, y_init = np.array(np.meshgrid(np.linspace(-0.01, 0.01, 21), np.linspace(-0.01, 0.01, 21))).reshape(2,-1)
    px_init = 0
    py_init = 0
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                   particle_ref=xc.fluka.engine.particle_ref,
                                   _capacity=xc.fluka.engine.capacity)

    part = part_init.copy()
    xc.fluka.engine.stop(clean=True)


def test_prototypes():
    raise ValueError("Need to write test for FlukaAssembly to check registry works as expected")

