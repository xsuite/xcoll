# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import pytest
import numpy as np
from pathlib import Path

import xtrack as xt
import xpart as xp
import xcoll as xc

from xobjects.test_helpers import for_all_test_contexts
try:
    import rpyc
except ImportError as e:
    rpyc = None


path = Path(__file__).parent / 'data'
particle_ref = xt.Particles('proton', p0c=6.8e12)


@pytest.mark.skipif(not xc.geant4.environment.ready, reason="BDSIM+Geant4 installation not found")
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
    _capacity = num_part*4
    coll = xc.Geant4Collimator(length=0.6, jaw=0.001, material='Ti')
    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.reentry_protection_enabled = False
    part = xp.build_particles(x=np.random.normal(coll.jaw_L + 1e-4, 1e-4, num_part),
                              particle_ref=particle_ref, _capacity=_capacity)

    cwd = Path.cwd()
    xc.geant4.engine.start(elements=coll, seed=1993)
    assert cwd != Path.cwd()
    temp_cwd = xc.geant4.engine._cwd
    assert temp_cwd == Path.cwd()
    assert xc.geant4.engine._already_started
    assert not Path('rpyc.log').exists()
    assert Path('root.out').exists()
    assert Path('root.err').exists()
    assert Path('geant4.out').exists()
    assert Path('geant4.err').exists()
    assert Path('engine.out').exists()
    assert Path('engine.err').exists()


    t_start = time.time()
    coll.track(part)
    print(f"Time per track: {(time.time()-t_start)*1e3:.2f}ms for "
        + f"{num_part} protons through {coll.length:.2f}m")
    assert (part.state == xc.headers.particle_states.LOST_WITHOUT_SPEC).sum() == 0   # No particles should be lost without specification
    assert (part.state == xc.headers.particle_states.LOST_ON_GEANT4_COLL).sum() > 0  # Some particles should have died in the collimator
    assert (part.state == 1).sum() > 0

    xc.geant4.engine.stop(clean=True)
    assert cwd == Path.cwd()
    assert not temp_cwd.exists()
    assert xc.geant4.engine._cwd is None
    assert xc.geant4.engine._g4link is None
    assert xc.geant4.engine._already_started
    assert not Path('rpyc.log').exists()
    assert not Path('root.out').exists()
    assert not Path('root.err').exists()
    assert not Path('geant4.out').exists()
    assert not Path('geant4.err').exists()
    assert not Path('engine.out').exists()
    assert not Path('engine.err').exists()

    with pytest.raises(RuntimeError, match="Cannot restart Geant4 engine in non-reentry-safe mode"):
        xc.geant4.engine.start(elements=coll, seed=1993)
