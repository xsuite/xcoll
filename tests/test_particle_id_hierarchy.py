# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import pytest
import numpy as np

import xpart as xp
import xtrack as xt
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
def test_geant4(test_context):
    length = 0.6
    num_slices = 1000
    coll = xc.Geant4Collimator(length=length/num_slices, jaw=0.001, material='Ti', _context=test_context)
    xc.geant4.engine.particle_ref = particle_ref
    # xc.geant4.engine.reentry_protection_enabled = False
    xc.geant4.engine.start(elements=coll, relative_energy_cut=0.1,
                           bdsim_config_file=path / f'geant4_protons.gmad')

    step = 2e-5
    _capacity = 25_000
    x = np.arange(0.002, 0.003+step, step)
    y = np.arange(-0.001, 0.001+step, step)
    X, Y = np.meshgrid(x,y)
    coords = np.vstack([X.ravel(), Y.ravel()]).T
    part_init = xp.build_particles(x=coords[:,0], y=coords[:,1], particle_ref=particle_ref,
                                   _capacity=_capacity, _context=test_context)
    parents = part_init.particle_id[part_init.state==1]
    # children = part_init.particle_id
    # children = children[children>=0]
    # children = {int(cc): [] for cc in children}

    part = part_init.copy()
    t_start = time.time()
    for _ in range(num_slices):
        coll.track(part)
    print(f"Time per track: {(time.time()-t_start)/num_slices*1e3:.2f}ms for "
        + f"{len(part_init.x)} protons through {coll.length:.2f}m")

    xc.geant4.engine.stop(clean=True)

    mask_child = part.particle_id > parents.max()
    assert not np.any(part.parent_particle_id[mask_child] == part.particle_id[mask_child])

    # mask = part.parent_particle_id != part.particle_id

    # xdiffs = []
    # ydiffs = []
    # pxvals = []
    # pyvals = []
    # for parent_id,x1,y1,px,py in zip(part.parent_particle_id[mask],part.x[mask],part.y[mask],part.px[mask],part.py[mask]):
    #     mask2 = part_init.particle_id == parent_id
    #     x0 = part_init.x[mask2][0]
    #     y0 = part_init.y[mask2][0]
    #     dx = x1 - x0
    #     dy = y1 - y0
    #     xdiffs.append(dx)
    #     ydiffs.append(dy)
    #     pxvals.append(px)
    #     pyvals.append(py)

    # assert np.max(np.abs(xdiffs)) < 5e-3
    # assert np.max(np.abs(ydiffs)) < 5e-3
