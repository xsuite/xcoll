# This is a script to test the preliminary FLUKA integration into xtrack.
# The FORTRAN source code is in xcoll/scattering_routines/fluka/FORTRAN_src/
#     (wrapper in pyfluka.f90, original SixTrack code in mod_fluka.f90 and others)
# The FlukaIO code is in xcoll/scattering_routines/fluka/flukaio/
#
# Before running the script, the source code needs to be compiled. In the root
# package folder, run ./compile_fluka.sh to do so.

import numpy as np
import xpart as xp
import xcoll as xc
import time


_capacity = 20000

coll = xc.FlukaCollimator(length=0.753)
coll.jaw = 0.001
coll.gap = 1  # STUB, needed to make FLUKAbuilder succeed but so far no clue why


# Connect to FLUKA
xc.FlukaEngine.start(elements=coll, names='tcp.c6l7.b1', debug_level=1, _capacity=_capacity)
particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
xc.FlukaEngine.set_particle_ref(particle_ref=particle_ref)

# Create an initial distribution of particles, random in 4D (with the
# longitudinal coordinates set to zero)
num_part = int(10000)
x_init   = np.random.normal(loc=1.288e-3, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init, particle_ref=particle_ref,
                          _capacity=_capacity)


# Do the tracking
start = time.time()
coll.track(part)
print(f"Tracking {num_part} particles took {round(time.time()-start,1)}s")


# Stop the FLUKA server
xc.FlukaEngine.stop()
