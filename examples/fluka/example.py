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


# Start rfluka and flukaserver and make the connection, based on
# the input files in this folder (lhc_run3_30cm.inp and insertion.txt)
xc.FlukaEngine.start_server("lhc_run3_30cm.inp")


# Create a FlukaCollimator beam element (ID 31 is the TCP.C6L7.B1 as
# defined by the input files)
coll = xc.FlukaCollimator(collimator_id=31, length=1.48200)


# Create an initial distribution of particles, random in 4D (with the
# longitudinal coordinates set to zero)
#num_part = int(10)
num_part = int(20000)
x_init   = np.random.normal(loc=1.288e-3, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
#part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)
part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, energy0=7e12)

# Do the tracking.
# This will:
#    1) prepare the arrays to be fortran-compatible            (OK)
#    2) call the fortran wrapper track_fluka() in pyfluka.f90  (OK)
#    3) which calls fluka_send_receive() in mod_fluka.f90      (OK)
#    4) which calls fluka_send(), then fluka_receive()
#npart, max_part, x_part, xp_part, y_part, yp_part, zeta_part, e_part, m_part, q_part, \
#A_part, Z_part, pdgid_part, partID, parentID, partWeight, spin_x_part, spin_y_part, spin_z_part = \
#    coll.track(part)
coll.track(part)

xc.FlukaEngine.stop_server()
