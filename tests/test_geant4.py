# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest
import numpy as np

import xtrack as xt
import xpart as xp
import xcoll as xc

from xobjects.test_helpers import for_all_test_contexts
try:
    import rpyc
except ImportError as e:
    rpyc = None

path = xc._pkg_root.parent / 'tests' / 'data'
particle_ref = xt.Particles('proton', p0c=6.8e12)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(not xc.geant4.environment.compiled, reason="BDSIM+Geant4 installation not found")
@pytest.mark.skipif(rpyc is None, reason="rpyc not installed")
def test_reload_bdsim(test_context):
    num_part = 1000
    coll = xc.Geant4Collimator(length=0.6, jaw=0.001, material='Ti', _context=test_context)
    xc.geant4.engine.particle_ref = particle_ref
    part_init = xp.build_particles(x=np.random.normal(coll.jaw_L, num_part),
                                   particle_ref=particle_ref, _capacity=4*num_part)
    part = []
    for _ in range(3):
        part.append(part_init.copy())
        xc.geant4.engine.start(elements=coll, seed=1993,
                               bdsim_config_file=path / 'geant4_protons.gmad')
        coll.track(part[-1])
        assert (part[-1].state == xc.headers.particle_states.LOST_WITHOUT_SPEC).sum() == 0   # No particles should be lost without specification
        assert (part[-1].state == xc.headers.particle_states.LOST_ON_GEANT4_COLL).sum() > 0  # Some particles should have died in the collimator
        assert (part[-1].state == 1).sum() > 0                                               # Some particles should have survived
        xc.geant4.engine.stop()

    # Check that the particles are the same
    for i in range(1, 3):
        assert np.allclose(part[0].x, part[i].x, atol=1e-12)
        assert np.allclose(part[0].px, part[i].px, atol=1e-12)
        assert np.allclose(part[0].y, part[i].y, atol=1e-12)
        assert np.allclose(part[0].py, part[i].py, atol=1e-12)
        assert np.allclose(part[0].zeta, part[i].zeta, atol=1e-12)
        assert np.allclose(part[0].delta, part[i].delta, atol=1e-12)
        assert np.all(part[0].state == part[i].state)
        assert np.all(part[0].particle_id == part[i].particle_id)
        assert np.all(part[0].parent_particle_id == part[i].parent_particle_id)


def test_serial_bdsim(pytestconfig):
    # Skip if Geant4Engine has already been started
    if xc.geant4.engine._already_started:
        pytest.skip("Cannot run serial BDSIM test - Geant4Engine has ran before.")

    # Verify that this test is not run in parallel
    n = None
    try:
        n = pytestconfig.getoption("numprocesses")  # xdist option
    except Exception:
        pass
    try:
        n = int(n)
    except Exception:
        if n:  # e.g. 'auto'
            n = 999999
    if n and n > 1:
        pytest.skip("Cannot run serial BDSIM test in parallel.")

    num_part = 1000
    coll = xc.Geant4Collimator(length=0.6, jaw=0.001, material='Ti')
    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.reentry_protection_enabled = False
    part = xp.build_particles(x=np.random.normal(coll.jaw_L, num_part),
                                   particle_ref=particle_ref, _capacity=4*num_part)
    xc.geant4.engine.start(elements=coll, seed=1993,
                           bdsim_config_file=path / 'geant4_protons.gmad')
    assert xc.geant4.engine._already_started
    coll.track(part)
    assert (part.state == xc.headers.particle_states.LOST_WITHOUT_SPEC).sum() == 0   # No particles should be lost without specification
    assert (part.state == xc.headers.particle_states.LOST_ON_GEANT4_COLL).sum() > 0  # Some particles should have died in the collimator
    assert (part.state == 1).sum() > 0
    xc.geant4.engine.stop()
    assert xc.geant4.engine._already_started
    with pytest.raises(RuntimeError, match="Cannot restart Geant4 engine in non-reentry-safe mode"):
        xc.geant4.engine.start(elements=coll, seed=1993,
                               bdsim_config_file=path / 'geant4_protons.gmad')


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(not xc.geant4.environment.compiled, reason="BDSIM+Geant4 installation not found")
@pytest.mark.skipif(rpyc is None, reason="rpyc not installed")
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
                                     _context=test_context, material='cu')
        g4_collimators.append(g4coll)
        bacoll = xc.BlackAbsorber(length=L, angle=angle, jaw=jaws+shift,
                                  _context=test_context)
        ba_collimators.append(bacoll)

    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=g4_collimators, seed=1993,
                           bdsim_config_file=path / f'geant4_ba_protons.gmad')

    x = np.random.uniform(-0.1, 0.1, n_part)
    y = np.random.uniform(-0.1, 0.1, n_part)
    px = np.random.uniform(-0.1, 0.1, n_part)
    py = np.random.uniform(-0.1, 0.1, n_part)
    part = xp.build_particles(x=x, y=y, px=px, py=py, _context=test_context,
                              particle_ref=xc.geant4.engine.particle_ref,
                              _capacity=2*n_part)
    part_ba = part.copy()
    # xc.geant4.engine.particle_ref

    for coll in g4_collimators:
        coll.track(part)

    for coll in ba_collimators:
        coll.track(part_ba)

    part.sort(interleave_lost_particles=True)
    part_ba.sort(interleave_lost_particles=True)

    assert np.all(part.filter(part.state==1).particle_id 
                  == part_ba.filter(part_ba.state==1).particle_id)

    # Stop the Geant4 connection
    xc.geant4.engine.stop(clean=True)
