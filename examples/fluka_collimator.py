# This is a script to test the FLUKA integration into xtrack.
# Before running the script, the source code needs to be compiled. In the root
# package folder, run ./compile_fluka.sh to do so.

import numpy as np
import xpart as xp
import xtrack as xt
import xcoll as xc
import time

import matplotlib.pyplot as plt


num_part = int(10_000)
_capacity = num_part*2


if xc.fluka.engine.is_running():
    xc.fluka.engine.stop(clean=True)


# Create a FLUKA collimator
coll = xc.FlukaCollimator(length=0.6, material='mogr')
coll.jaw = 0.001

# The same collimator in Everest
coll2 = xc.EverestCollimator(length=0.6, material=xc.materials.MolybdenumGraphite)
coll2.jaw = 0.001


# Connect to FLUKA
xc.fluka.engine.particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12)
xc.fluka.engine.capacity = _capacity
xc.fluka.engine.seed = 5656565
xc.fluka.engine.start(elements=coll, clean=True, verbose=False)


# Create an initial distribution of particles, random in 4D, on the left jaw (with the
# longitudinal coordinates set to zero)
np.random.seed(seed=23823842)
x_init   = np.random.normal(loc=0.001, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                               particle_ref=xc.fluka.engine.particle_ref,
                               _capacity=xc.fluka.engine.capacity)
part = part_init.copy()
part2 = part_init.copy()
part_test = part_init.copy()


# Do the tracking in FLUKA
# coll.track(part_test)  # pre-track to compile the code for a fair comparison
print(f"Tracking {num_part} particles (FLUKA)...     ", end='')
start = time.time()
coll.track(part)
print(f"Done in {round(time.time()-start, 3)}s.")

# Do the tracking in Everest
# coll2.track(part_test)  # pre-track to compile the code for a fair comparison
print(f"Tracking {num_part} particles (Everest)...   ", end='')
start = time.time()
coll2.track(part2)
print(f"Done in {round(1000*(time.time()-start), 3)}ms")

print(f"Survived in FLUKA: {len(part.state[part.state>0])}/{num_part}")
print(f"Survived in Everest: {len(part2.state[part2.state>0])}/{num_part}")


# Stop the FLUKA server
xc.fluka.engine.stop(clean=True)


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
