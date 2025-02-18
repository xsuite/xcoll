import numpy as np
import pytest
import time
from scipy.stats import ks_2samp

import xobjects as xo
import xpart as xp
import xcoll as xc

from xobjects.test_helpers import for_all_test_contexts

# try the import here and skip tests if missing
# also need the import here in case of pytest --forked
try:
    import collimasim as cs
except ImportError:
    cs = None

path = xc._pkg_root.parent / 'tests' / 'data'


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_black_absorbers(test_context):
    n_part = 10000
    angles=[0,45,90]
    jaws = [0.03, -0.02]
    co = [-0.01, 0.01]
    L = 0.873

    g4_collimators = []
    ba_collimators = []
    for ii, angle in enumerate(angles):
        shift = co[0]*np.cos(angle) + co[1]*np.sin(angle)
        g4coll = xc.Geant4Collimator(length=L, angle=angle, jaw=jaws+shift,
                                     _context=test_context, material='cu',
                                     geant4_id=f'g4coll_{ii}')
        g4_collimators.append(g4coll)
        bacoll = xc.BlackAbsorber(length=L, angle=angle, jaw=jaws+shift,
                                  _context=test_context)
        ba_collimators.append(bacoll)

    xc.Geant4Engine.start(elements=g4_collimators, seed=1993, particle_ref='proton', p0c=7e12,
                          bdsim_config_file=str(path / f'geant4_ba_protons.gmad'))

    x = np.random.uniform(-0.1, 0.1, n_part)
    y = np.random.uniform(-0.1, 0.1, n_part)
    px = np.random.uniform(-0.1, 0.1, n_part)
    py = np.random.uniform(-0.1, 0.1, n_part)
    part = xp.build_particles(x=x, y=y, px=px, py=py, _context=_context,
                              particle_ref=xc.Geant4Engine().particle_ref)
    part_ba = part.copy()
    xc.Geant4Engine().particle_ref

    for coll in g4_collimators:
        coll.track(part)

    for coll in ba_collimators:
        coll.track(part_ba)

    part.sort(interleave_lost_particles=True)
    part_ba.sort(interleave_lost_particles=True)

    assert np.all(part.filter(part.state==1).particle_id 
                  == part_ba.filter(part_ba.state==1).particle_id)


@pytest.mark.parametrize('num_part', [1000, 5000])
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_simple_track(num_part):
    print(f"Running test_simple_track with {num_part} particles")
    _capacity = num_part*2

    # If a previous test failed, stop the server manually
    if xc.Geant4Engine.is_running():
        xc.Geant4Engine.stop(clean=True)

    # Define Geant4 collimator and start Geant4 engine
    coll = xc.Geant4Collimator(material='mogr', length=0.6, jaw=0.001)
    xc.Geant4Engine.start(elements=coll, seed=1993,
                          bdsim_config_file=path / 'geant4_protons.gmad')

    # Particle distribution
    x_init   = np.random.normal(loc=0.002, scale=1e-3, size=num_part)
    px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                   particle_ref=particle_ref, _capacity=_capacity)
    part_geant4 = part_init.copy()
    part_drift = part_init.copy()

    # Geant4 tracking
    start = time.time()
    coll.track(part_geant4)
    geant4_time = round(time.time()-start, 3)
    expected_time = round(0.025 * num_part)    # To be adapted to reasonable Geant4 time
    if geant4_time > expected_time:
        raise RuntimeError(f"Geant4 tracking took {geant4_time}s, exceeding threshold of {expected_time}s")
    else:
        print(f"Geant4 tracking took {geant4_time}s")

    # Drift tracking
    coll._equivalent_drift.length = coll.length
    coll._equivalent_drift.track(part_drift)

    mask_geant4 = part_geant4.state > 0
    mask_drift = part_drift.state > 0
    ks_stat, p_value = ks_2samp(part_geant4.x[mask_geant4], part_drift.x[mask_drift])
    assert p_value <= 0.05, f"Distributions are not significantly different (p = {p_value})"
    print(f"KS test passed with p = {p_value}")
    mask_geant4 = part_geant4.state > 0
    mask_drift = part_drift.state > 0
    ks_stat, p_value = ks_2samp(part_geant4.x[mask_geant4], part_drift.x[mask_drift])
    assert p_value <= 0.05, f"Distributions are not significantly different (p = {p_value})"
    print(f"KS test passed with p = {p_value}")

    # Stop the Geant4 connection
    xc.Geant4Engine.stop(clean=True)


# @pytest.mark.parametrize('tilt', [0, [2.2e-6, 1.3e-6], [1.9e-6, -2.7e-6]],
#                          ids=['no_tilt', 'positive_tilt', 'pos_neg_tilt'])
@pytest.mark.parametrize('angle', [0, 90, 130.5])
@pytest.mark.parametrize('jaw', [0.001, [0.0013, -0.002789], [-1.2e-6, -3.2e-3], [3.789e-3, 4.678e-7]],
                         ids=['symmetric', 'asymmetric', 'negative', 'positive'])
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_jaw(jaw, angle): #, tilt):
    tilt = 0    # For now, need to implement test for tilted jaws
    _ACCURACY = 1e-12  # Anything in this region around the jaw might or might not hit; we can't be sure
    num_part = 5000
    _capacity = num_part*2
    jaw_band = 1.e-6

    # If a previous test failed, stop the server manually
    if xc.Geant4Engine.is_running():
        xc.Geant4Engine.stop(clean=True)

    # Define Geant4 collimator and start Geant4 engine
    coll = xc.Geant4Collimator(material='mogr', length=0.6, jaw=jaw, angle=angle, tilt=tilt)
    xc.Geant4Engine.start(elements=coll, seed=1993,
                          bdsim_config_file=path / 'geant4_protons.gmad')
    particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)

    # Particle distribution (x and y are in the frame of the collimator)
    num_part_step = num_part//5
    x = np.random.uniform(-0.02, 0.02, num_part_step)
    x = np.concatenate([x, np.random.uniform(coll.jaw_L - jaw_band, coll.jaw_L -_ACCURACY, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_L +_ACCURACY, coll.jaw_L + jaw_band, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_R - jaw_band, coll.jaw_R -_ACCURACY, num_part_step)])
    x = np.concatenate([x, np.random.uniform(coll.jaw_R +_ACCURACY, coll.jaw_R + jaw_band, num_part_step)])
    y = np.random.uniform(-0.02, 0.02, 5*num_part_step)
    x_new = np.cos(np.deg2rad(angle))*x - np.sin(np.deg2rad(angle))*y
    y_new = np.sin(np.deg2rad(angle))*x + np.cos(np.deg2rad(angle))*y
    part_init = xp.build_particles(x=x_new, y=y_new, particle_ref=xc.Geant4Engine().particle_ref,
                                   _capacity=xc.Geant4Engine()._capacity)
    mask = np.concatenate([(x >= min(coll.jaw_LU, coll.jaw_LD)) | (x <= max(coll.jaw_RU, coll.jaw_RD)),
                          np.full(5*num_part_step, False)])
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
    # plt.axvline(coll.jaw_L, c='k', linestyle='--')
    # plt.axvline(coll.jaw_R, c='k', linestyle='--')
    # plt.show()


    # Particles that are supposed to not have hit the collimator, but have a kick or are dead, are considered faulty
    assert not np.any(abs(part.px[mask_not_hit]) > _ACCURACY)
    assert not np.any(abs(part.py[mask_not_hit]) > _ACCURACY)
    assert not np.any(part.state[mask_not_hit] < 1)

    # Particles that are supposed to have hit the collimator, but are alive and have no kick, are considered faulty
    faulty =  mask_hit & (abs(part.px) < _ACCURACY) & (abs(part.py) < _ACCURACY)
    faulty &= (part.state > 0)
    assert len(part.x[faulty]) <= 1  # We allow for a small margin of error

    # Stop the Geant4 connection
    xc.Geant4Engine.stop(clean=True)
