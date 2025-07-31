# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import json
import numpy as np
from pathlib import Path
import xpart as xp
import xtrack as xt
import xcoll as xc
import xobjects as xo
import pytest
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts
import sys
import time

# try the import here and skip tests if missing
# also need the import here in case of pytest --forked
try:
    import collimasim as cs
except ImportError:
    cs = None


path = xc._pkg_root.parent / 'tests' / 'data'
particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
@pytest.mark.skipif(cs is None, reason="Geant4 tests need collimasim installed")
def test_geant4_ppid(test_context):
    coll = xc.Geant4Collimator(length=0.6, jaw=0.001, material='Ti')
    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=coll, relative_energy_cut=0.1,
                          bdsim_config_file=path / f'geant4_protons.gmad')

    num_part_1d = 1000
    jaw_band = 1e-6
    _capacity = int(2e6)
    x = np.linspace(0.002, 10, num_part_1d)
    y = np.linspace(-10, 10, num_part_1d)
    X, Y = np.meshgrid(x,y)
    coords = np.vstack([X.ravel(), Y.ravel()]).T
    part_init = xp.build_particles(x=coords[:,0], y=coords[:,1], particle_ref=particle_ref, _capacity=_capacity)

    part = part_init.copy()
    coll.track(part)
    xc.geant4.engine.stop()
 
    mask = part.parent_particle_id != part.particle_id

    xdiffs = []
    ydiffs = []
    pxvals = []
    pyvals = []
    for parent_id,x1,y1,px,py in zip(part.parent_particle_id[mask],part.x[mask],part.y[mask],part.px[mask],part.py[mask]):
        mask2 = part_init.particle_id == parent_id
        x0 = part_init.x[mask2][0]
        y0 = part_init.y[mask2][0]
        dx = x1 - x0
        dy = y1 - y0
        xdiffs.append(dx)
        ydiffs.append(dy)
        pxvals.append(px)
        pyvals.append(py)

    assert np.max(np.abs(xdiffs)) < 5e-3
    assert np.max(np.abs(ydiffs)) < 5e-3
