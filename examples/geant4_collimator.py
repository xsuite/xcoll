# This is a script to test the Geant4 integration into xtrack.
# Before running the script, the source code needs to be compiled.
# See the geant4_init.py example to do so.

import numpy as np
import xpart as xp
import xtrack as xt
import xcoll as xc
import time

import matplotlib.pyplot as plt


num_part = int(10_000)
_capacity = num_part*2


if xc.geant4.engine.is_running():
    xc.geant4.engine.stop(clean=True)


# Create a Geant4 collimator
coll = xc.Geant4Collimator(length=0.6, material='mogr')
coll.jaw = 0.001

# The same collimator in Everest
coll2 = xc.EverestCollimator(length=0.6, material=xc.materials.MolybdenumGraphite)
coll2.jaw = 0.001


# Connect to Geant4
xc.geant4.engine.particle_ref = xt.Particles('proton', p0c=6.8e12)
xc.geant4.engine.seed = 5656565
xc.geant4.engine.start(elements=coll, clean=True, verbose=False)


# Create an initial distribution of particles, random in 4D, on the left jaw (with the
# longitudinal coordinates set to zero)
np.random.seed(seed=23823842)
x_init   = np.random.normal(loc=0.001, scale=0.2e-3, size=num_part)
px_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                               particle_ref=xc.geant4.engine.particle_ref,
                               _capacity=_capacity)
part = part_init.copy()
part2 = part_init.copy()
part_test = part_init.copy()


# Do the tracking in Geant4
coll.track(part_test)  # pre-track to compile the code for a fair comparison
print(f"Tracking {num_part} particles (Geant4)...     ", end='')
start = time.time()
coll.track(part)
print(f"Done in {round(time.time()-start, 3)}s.")

# Do the tracking in Everest
coll2.track(part_test)  # pre-track to compile the code for a fair comparison
print(f"Tracking {num_part} particles (Everest)...   ", end='')
start = time.time()
coll2.track(part2)
print(f"Done in {round(1000*(time.time()-start), 3)}ms")

print(f"Survived in Geant4: {len(part.state[part.state>0])}/{num_part}")
print(f"Survived in Everest: {len(part2.state[part2.state>0])}/{num_part}")


# Stop the Geant4 server
xc.geant4.engine.stop(clean=True)


# Plot distribution
cutoff = 0.0001
_, ax = plt.subplots(2, 2, figsize=(12, 9))
mask = part.state > 0
parents = part.parent_particle_id[mask]
output = part.kin_xprime[mask] - part_init.kin_xprime[parents]
new_mask = (output > -cutoff) & (output < cutoff)   # For visibility
ax[0, 0].hist2d(part_init.x[parents][new_mask], output[new_mask], 200)
ax[0, 0].axhline(coll.jaw_L, color='black', ls='--')
# ax[0, 0].set_xlabel('x [m]')
ax[0, 0].set_ylabel(r'$\theta_f$ - $\theta_i$ [rad]')
ax[0, 0].set_ylim(-cutoff, cutoff)
ax[0, 0].set_title('Geant4')

mask2 = part2.state > 0
parents2 = part2.parent_particle_id[mask2]
output2 = part2.kin_xprime[mask2] - part_init.kin_xprime[parents2]
new_mask2 = (output2 > -cutoff) & (output2 < cutoff)
ax[0, 1].hist2d(part_init.x[parents2][new_mask2], output2[new_mask2], 200)
ax[0, 1].axhline(coll2.jaw_L, color='black', ls='--')
ax[0, 1].set_ylim(-cutoff, cutoff)
ax[0, 1].set_title('Everest')

# Plot only particles that hit
new_mask = new_mask & (part_init.x[parents] >= coll.jaw_L)
ax[1, 0].hist2d(part_init.x[parents][new_mask], output[new_mask], 200)
ax[1, 0].axhline(coll.jaw_L, color='black', ls='--')
ax[0, 1].set_ylim(-cutoff, cutoff)
ax[1, 0].set_xlabel('x [m]')
ax[1, 0].set_ylabel(r'$\theta_f$ - $\theta_i$ [rad]')

new_mask2 = new_mask2 & (part_init.x[parents2] >= coll2.jaw_L)
ax[1, 1].hist2d(part_init.x[parents2][new_mask2], output2[new_mask2], 200)
ax[1, 1].axhline(coll2.jaw_L, color='black', ls='--')
ax[0, 1].set_ylim(-cutoff, cutoff)
ax[1, 1].set_xlabel('x [m]')

plt.tight_layout()
plt.show()
