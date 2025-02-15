import numpy as np
import pytest
import time
from scipy.stats import ks_2samp

import xobjects as xo
import xpart as xp
import xcoll as xc

from xobjects.test_helpers import for_all_test_contexts

#import pydevd
#pydevd.settrace(suspend=True)
#pydevd.settrace('127.0.0.1', port=5678, stdoutToServer=True, stderrToServer=True, suspend=True)

# try the import here and skip tests if missing
# also need the import here in case of pytest --forked
try:
    import collimasim as cs
except ImportError:
    cs = None

path = xc._pkg_root.parent / 'tests' / 'data'

jaw = 0.001
#jaw = [0.0013, -0.002789]
#jaw = [-1.2e-6, -3.2e-3]
angle = 0

tilt = 0    # For now, need to implement test for tilted jaws
_ACCURACY = 1.5e-12  # Anything in this region around the jaw might or might not hit; we can't be sure
num_part = 50
_capacity = num_part*2
jaw_band = 4.e-9
geant4_margin = 2.e-9

# If a previous test failed, stop the server manually
if xc.Geant4Engine.is_running():
    xc.Geant4Engine.stop(clean=True)

# Define collimator and start the Geant4 server
coll = xc.Geant4Collimator(length=0.6, jaw=0.001, angle=0, tilt=0, material='Ti', geant4_id=f'g4coll_0')


# Particle distribution (x and y are in the frame of the collimator)
num_part_step = num_part//5
context = xo.ContextCpu()
particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=7e12, _context=context)
x = np.random.uniform(-0.02, 0.02, num_part_step)
x = np.concatenate([x, np.random.uniform(coll.jaw_L - jaw_band, coll.jaw_L -geant4_margin -_ACCURACY, num_part_step)])
x = np.concatenate([x, np.random.uniform(coll.jaw_L +geant4_margin +_ACCURACY, coll.jaw_L + jaw_band, num_part_step)])
x = np.concatenate([x, np.random.uniform(coll.jaw_R - jaw_band, coll.jaw_R -geant4_margin -_ACCURACY, num_part_step)])
x = np.concatenate([x, np.random.uniform(coll.jaw_R +geant4_margin +_ACCURACY, coll.jaw_R + jaw_band, num_part_step)])
y = np.random.uniform(-0.02, 0.02, 5*num_part_step)
x_new = np.cos(np.deg2rad(angle))*x - np.sin(np.deg2rad(angle))*y
y_new = np.sin(np.deg2rad(angle))*x + np.cos(np.deg2rad(angle))*y
part_init = xp.build_particles(x=x_new, y=y_new, particle_ref=xc.Geant4Engine().particle_ref,
                               _capacity=_capacity)

mask = np.concatenate([(x >= min(coll.jaw_LU, coll.jaw_LD)) | (x <= max(coll.jaw_RU, coll.jaw_RD)),
                      np.full(5*num_part_step, False)])
hit_ids = part_init.particle_id[mask & (part_init.state > 0)]
not_hit_ids = part_init.particle_id[~mask & (part_init.state > 0)]

# TODO: jaw tilts, and particle angles

# Track
part = part_init.copy()
if True:
    xc.Geant4Engine.start(elements=coll, seed=1993, particle_ref='proton', p0c=7.e12,
                          bdsim_config_file=str(path / 'geant4_protons.gmad'))
    coll.track(part)
    xc.Geant4Engine.stop()

    xc.Geant4Engine.start(elements=coll, seed=1993, particle_ref='proton', p0c=7.e12,
                          bdsim_config_file=str(path / 'geant4_protons.gmad'))
    part2 = part_init.copy()
    coll.track(part2)
    xc.Geant4Engine.stop()

xc.Geant4Engine.start(elements=coll, seed=1992, particle_ref='proton', p0c=7.e12,
                      bdsim_config_file=str(path / 'geant4_protons.gmad'))
part3 = part_init.copy()
coll.track(part3)
print(part_init.px)
print(part3.px)
print()
print(part_init.state)
print(part3.state)
print()
print(part_init.parent_particle_id)
print(part3.parent_particle_id)
print('hej3')

