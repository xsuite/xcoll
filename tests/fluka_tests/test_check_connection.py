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

@pytest.fixture
def fluka_collimator(coll_jaw):
    # length = 0.6 # XXX TO BE REMOVED
    coll = xc.FlukaCollimator(length=0.6)
    # coll_name = 'tcp.c6l7.b1'
    coll.jaw = coll_jaw
    # coll.gap = 1 # XXX STILL NEEDED?
    return coll

@pytest.fixture
def particle_distribution(num_part):
    x_init   = np.random.normal(loc=0.002, scale=1e-3, size=num_part)
    px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    xc.FlukaEngine.set_particle_ref(particle_ref=particle_ref)
    return xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init, particle_ref=particle_ref, _capacity=num_part*2)

def perform_ks_test(part1, part2):
    mask_fluka = part1.state > 0
    mask_drift = part2.state > 0
    ks_stat, p_value = ks_2samp(part1.x[mask_fluka], part2.x[mask_drift])
    assert p_value <= 0.05, f"Distributions are not significantly different (p = {p_value})"
    print(f"KS test passed with p = {p_value}")

@pytest.mark.parametrize("num_part", [1000])
@pytest.mark.parametrize("fluka_time_threshold, drift_time_threshold", [(25, 1)])
@pytest.mark.parametrize("coll_jaw", [0.001])
def test_fluka_connection(fluka_collimator, particle_distribution, num_part, fluka_time_threshold, drift_time_threshold, coll_jaw):
    _capacity = num_part*2
    coll = fluka_collimator

    # Connect to FLUKA
    coll_name = 'tcp.c6l7.b1'
    xc.FlukaEngine.start(elements=coll, names=coll_name, debug_level=1, _capacity=_capacity)

    # Particle distribution
    part_init = particle_distribution
    part_fluka = part_init.copy()
    part_drift = part_init.copy()

    # FLUKA tracking
    start = time.time()
    coll.track(part_fluka)
    fluka_time = round(time.time()-start, 3)
    assert fluka_time < fluka_time_threshold, f"FLUKA tracking took {fluka_time}s, exceeding threshold of {fluka_time_threshold}s with coll_jaw = {coll_jaw}" 
    
    # Drift tracking
    start = time.time()
    coll._equivalent_drift.length = 0.6 # length of tcp
    coll._equivalent_drift.track(part_drift)
    drift_time = round(time.time()-start, 3)
    assert drift_time < drift_time_threshold, f"Drift tracking took {drift_time}s, exceeding threshold of {drift_time_threshold}s with coll_jaw = {coll_jaw}"

    # Kolmogorov-Smirnov test to verify distribution difference

    # ks_stat, p_value = ks_2samp(part_fluka.x, part_drift.x)

    # assert p_value <= 0.05, f"Distributions are not significantly different (p = {p_value})"

    perform_ks_test(part_fluka, part_drift)

    # Stop the FLUKA server
    xc.FlukaEngine.stop()

    # plt.hist(part.x[mask], bins=50)
    # plt.savefig('fluka_hist.png')
    # plt.clf()
    # plt.scatter(part.x[mask], part.y[mask], label='FLUKA')
    # plt.savefig('fluka.png')
    # plt.clf()
    # plt.scatter(part2.x[mask2], part2.y[mask2], label='Drift')
    # plt.savefig('drift.png')
    # plt.show()


#     print(f"Survived in FLUKA: {len(part.state[part