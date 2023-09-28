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


# Start rfluka and flukaserver and make the connection, based on
# the input files in this folder (lhc_run3_30cm.inp and insertion.txt)
# The value in insertion.txt is the value to be used for inactive_front/back each
xc.FlukaEngine.start_server("lhc_run3_30cm.inp", fluka_ids={'tcp.c6l7.b1': 31}, n_alloc=21000)


# Create a FlukaCollimator beam element (ID 31 is the TCP.C6L7.B1 as
# defined by the input files; the length is the active_length + the value in insertion.txt)
coll = xc.FlukaCollimator(fluka_id=31, length=1.48200)


# Set a reference particle
particle_ref = xp.Particles.build_reference_particle(pdg_id='proton', p0c=7e12)
xc.FlukaEngine().set_particle_ref(particle_ref)


# Create an initial distribution of particles, random in 4D (with the
# longitudinal coordinates set to zero)
num_part = int(10000)
x_init   = np.random.normal(loc=1.288e-3, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init, particle_ref=particle_ref, _capacity=20000)


# Do the tracking
start = time.time()
coll.track(part)
print(f"Tracking {num_part} particles took {round(time.time()-start,1)}s")


# Stop the FLUKA server
xc.FlukaEngine.stop_server()
