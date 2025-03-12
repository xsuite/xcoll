import numpy as np
import xpart as xp
import xcoll as xc
import time

import matplotlib.pyplot as plt

num_part = int(100)
_capacity = num_part*1000

if xc.FlukaEngine.is_running():
    xc.FlukaEngine.stop()

# Create a FLUKA collimator
coll = xc.FlukaCollimator(length=0.4, assembly='fcc_tcp')
coll.jaw = 0.001


# Connect to FLUKA
xc.FlukaEngine.particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='electron', p0c=200e9)
xc.FlukaEngine.capacity = _capacity
xc.FlukaEngine.start(elements=coll, clean=False, verbose=True, return_all=True,
                     return_neutral=True, electron_lower_momentum_cut=1.e9)


# Create an initial distribution of particles, random in 4D, on the left jaw (with the
# longitudinal coordinates set to zero)
x_init   = np.random.normal(loc=0.0015, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                               particle_ref=xc.FlukaEngine.particle_ref,
                               _capacity=xc.FlukaEngine.capacity)
part = part_init.copy()


# Do the tracking in FLUKA
print(f"Tracking {num_part} particles (FLUKA)...     ", end='')
start = time.time()
coll.track(part)
print(f"Done in {round(time.time()-start, 3)}s.")
print(f"Survived in FLUKA: {len(part.state[part.state>0])}/{num_part}")

# Stop the FLUKA server
xc.FlukaEngine.stop(clean=False)
