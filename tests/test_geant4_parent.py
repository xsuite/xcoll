# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
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
g4coll = xc.Geant4Collimator(length=0.6, jaw=0.001, angle=0, tilt=0,
                            material='Ti', geant4_id=f'g4coll_0')
xc.Geant4Engine.start(elements=g4coll, seed=1993, particle_ref='proton', p0c=7e12, relative_energy_cut=0.1,
                          bdsim_config_file=str(path / f'geant4_protons.gmad'))

num_part_step = int(1e3)
jaw_band = 1e-6
angle = 0
_capacity = int(2e6)
particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=7e12)
x = np.linspace(0.002,10,num_part_step)
y = np.linspace(-10,10,num_part_step)
X, Y = np.meshgrid(x,y)
coords = np.vstack([X.ravel(), Y.ravel()]).T
part_init = xp.build_particles(x=coords[:,0], y=coords[:,1], particle_ref=particle_ref, _capacity=_capacity)

part = part_init.copy()
t00 = time.time()
g4coll.track(part)
t01 = time.time() 

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