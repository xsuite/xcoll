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

path = Path('/home/a20/ドキュメント/work/git/xtrack/dev/xcoll_geant4/xcoll/tests/data')

g4coll = xc.Geant4Collimator(length=0.6, jaw=0.001, angle=0, tilt=0,
                            material='Ti', geant4_id=f'g4coll_0')
xc.Geant4Engine.start(elements=g4coll, seed=1993, particle_ref='proton', p0c=7e12, relative_energy_cut=0.1,
                          bdsim_config_file=str(path / f'geant4_protons.gmad'))

num_part_step = int(25e3)
jaw_band = 1e-6
angle = 0
_capacity = int(100e3)
particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=7e12)
x = np.random.uniform(g4coll.jaw_L + jaw_band, g4coll.jaw_L + 10*jaw_band, num_part_step)
x = np.concatenate([x, np.random.uniform(g4coll.jaw_R - jaw_band, g4coll.jaw_R - 10*jaw_band, num_part_step)])
y = np.random.uniform(-0.02, 0.02, 2*num_part_step)
x_new = np.cos(np.deg2rad(angle))*x - np.sin(np.deg2rad(angle))*y
y_new = np.sin(np.deg2rad(angle))*x + np.cos(np.deg2rad(angle))*y
part_init = xp.build_particles(x=x_new, y=y_new, particle_ref=particle_ref, _capacity=_capacity)

# Track
part = part_init.copy()
t00 = time.time()
g4coll.track(part)
t01 = time.time()

# Particle distribution (x and y are in the frame of the collimator)
num_part_step = int(25e3)
jaw_band = 1e-6
angle = 0
_capacity = int(100e3)
particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=7e12)
x = np.random.uniform(g4coll.jaw_L - jaw_band, g4coll.jaw_L - 10*jaw_band, num_part_step)
x = np.concatenate([x, np.random.uniform(g4coll.jaw_R + jaw_band, g4coll.jaw_R + 10*jaw_band, num_part_step)])
y = np.random.uniform(-0.02, 0.02, 2*num_part_step)
x_new = np.cos(np.deg2rad(angle))*x - np.sin(np.deg2rad(angle))*y
y_new = np.sin(np.deg2rad(angle))*x + np.cos(np.deg2rad(angle))*y
part_init1 = xp.build_particles(x=x_new, y=y_new, particle_ref=particle_ref, _capacity=_capacity)

# Track
part1 = part_init1.copy()
t10 = time.time()
g4coll.track(part1)
t11 = time.time()

assert np.sum(part.state == -333) > 40000 ## some particles should have died
assert np.sum(np.abs((part.px[part_init.state==1] - part_init.px[part_init.state==1])) == 0) == 0 ## all particles should have some kick in x and y
assert np.sum(np.abs((part.py[part_init.state==1] - part_init.py[part_init.state==1])) == 0) == 0 ## all particles should have some kick in x and y

assert np.sum((part1.state < 1) & (part1.state > -9999)) == 0 ## no particles should have died
assert np.sum(np.abs((part1.px[part_init1.state==1] - part_init1.px[part_init1.state==1])) > 0) == 0 ## no particles should have kicks in x and y
assert np.sum(np.abs((part1.py[part_init1.state==1] - part_init1.py[part_init1.state==1])) > 0) == 0 ## no particles should have kicks in x and y
