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

import matplotlib.pyplot as plt

num_part = int(10_000)
_capacity = num_part*2

if xc.FlukaEngine.is_running():
    xc.FlukaEngine.stop()

# Create a FLUKA collimator
coll = xc.FlukaCollimator(length=0.6, assembly='lhc_tcp')
coll_name = 'tcp.c6l7.b1'
coll.jaw = 0.001

# The same collimator in Everest
coll2 = xc.EverestCollimator(length=0.6, material=xc.materials.MolybdenumGraphite)
coll2.jaw = 0.001


# Connect to FLUKA
xc.FlukaEngine.particle_ref = xp.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
xc.FlukaEngine.capacity = _capacity
xc.FlukaEngine.start(elements=coll, names=coll_name, clean=True, verbose=False)


# Create an initial distribution of particles, random in 4D, on the left jaw (with the
# longitudinal coordinates set to zero)
x_init   = np.random.normal(loc=0.001, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                               particle_ref=xc.FlukaEngine.particle_ref,
                               _capacity=xc.FlukaEngine.capacity)
part = part_init.copy()
part2 = part_init.copy()


# Do the tracking
start = time.time()
coll.track(part)
print(f"Tracking {num_part} particles (FLUKA) took {round(time.time()-start, 3)}s")
start = time.time()
coll2.track(part2)
print(f"Tracking {num_part} particles (Everest) took {round(time.time()-start, 3)}s")

print(f"Survived in FLUKA: {len(part.state[part.state>0])}/{num_part}")
print(f"Survived in Everest: {len(part2.state[part2.state>0])}/{num_part}")


# Stop the FLUKA server
xc.FlukaEngine.stop(clean=True)


# Make some plots

# mask = part.state > -9999999
mask = part.state > 0
parents = part.parent_particle_id[mask]
# print(sorted(np.unique((part.x[mask] - part_init.x[parents])/part_init.kin_xprime[parents])))

output = part.kin_xprime[mask] - part_init.kin_xprime[parents]
new_mask = (output > -0.0001) & (output < 0.0001)
_ = plt.hist2d(part_init.x[parents][new_mask], output[new_mask], 200)
plt.show()


# mask = part.state > -9999999
mask2 = part2.state > 0
parents2 = part2.parent_particle_id[mask2]
# print(sorted(np.unique((part.x[mask] - part_init.x[parents])/part_init.kin_xprime[parents])))

output2 = part2.kin_xprime[mask2] - part_init.kin_xprime[parents2]
new_mask2 = (output2 > -0.0001) & (output2 < 0.0001)
_ = plt.hist2d(part_init.x[parents2][new_mask2], output2[new_mask2], 200)
plt.show()


new_mask = new_mask & (part_init.x[parents] >= coll.jaw_L)
_ = plt.hist2d(part_init.x[parents][new_mask], output[new_mask], 200)
plt.show()


new_mask2 = new_mask2 & (part_init.x[parents2] >= coll2.jaw_L)
_ = plt.hist2d(part_init.x[parents2][new_mask2], output2[new_mask2], 200)
plt.show()
