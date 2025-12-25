# This is a script to test the FLUKA integration into xtrack.
# Before running the script, the source code needs to be compiled. In the root
# package folder, run ./compile_fluka.sh to do so.

import numpy as np
import xpart as xp
import xtrack as xt
import xcoll as xc
import time

import matplotlib.pyplot as plt


num_part = 10_000
capacity = 2*num_part
particle_ref = xt.Particles('proton', p0c=6.8e12)

if xc.fluka.engine.is_running():
    xc.fluka.engine.stop(clean=True)

# Create a FLUKA collimator
coll1 = xc.FlukaCollimator(length=0.6, material='mogr', jaw=0.001)

# The same collimator in Everest
coll2 = xc.EverestCollimator(length=0.6, material=xc.materials.MolybdenumGraphite, jaw=0.001)


# Connect to FLUKA
xc.fluka.engine.particle_ref = particle_ref
xc.fluka.engine.capacity = capacity
xc.fluka.engine.seed = 5656565
xc.fluka.engine.start(elements=coll1, clean=True, verbose=False, fortran_debug_level=1)


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
part1 = part_init.copy()
part2 = part_init.copy()
part_pre = xp.build_particles(x=[0], particle_ref=xc.fluka.engine.particle_ref, _capacity=2)


# Do the tracking in FLUKA
print(f"Tracking {num_part} particles (FLUKA)...     ", end='')
start = time.time()
coll1.track(part1)
print(f"Done in {round(time.time()-start, 3)}s.")

# Do the tracking in Everest
print(f"Tracking {num_part} particles (Everest)...   ", end='')
start = time.time()
coll2.track(part2)
print(f"Done in {round(1000*(time.time()-start), 3)}ms")

print(f"Survived in FLUKA: {len(part1.state[part1.state>0])}/{num_part}")
print(f"Survived in Everest: {len(part2.state[part2.state>0])}/{num_part}")

# Stop the FLUKA server
xc.fluka.engine.stop(clean=True)


# Plot distribution
cutoff = 0.0001
_, ax = plt.subplots(2, 2, figsize=(12, 9))
mask = part1.state > 0
parents = part1.parent_particle_id[mask]
output = part1.kin_xprime[mask] - part_init.kin_xprime[parents]
new_mask = (output > -cutoff) & (output < cutoff)   # For visibility
ax[0, 0].hist2d(part_init.x[parents][new_mask], output[new_mask], 200)
ax[0, 0].axhline(coll1.jaw_L, color='black', ls='--')
# ax[0, 0].set_xlabel('x [m]')
ax[0, 0].set_ylabel(r'$\theta_f$ - $\theta_i$ [rad]')
ax[0, 0].set_ylim(-cutoff, cutoff)
ax[0, 0].set_title('FLUKA')

mask2 = part2.state > 0
parents2 = part2.parent_particle_id[mask2]
output2 = part2.kin_xprime[mask2] - part_init.kin_xprime[parents2]
new_mask2 = (output2 > -cutoff) & (output2 < cutoff)
ax[0, 1].hist2d(part_init.x[parents2][new_mask2], output2[new_mask2], 200)
ax[0, 1].axhline(coll2.jaw_L, color='black', ls='--')
ax[0, 1].set_ylim(-cutoff, cutoff)
ax[0, 1].set_title('Everest')

# Plot only particles that hit
new_mask = new_mask & (part_init.x[parents] >= coll1.jaw_L)
ax[1, 0].hist2d(part_init.x[parents][new_mask], output[new_mask], 200)
ax[1, 0].axhline(coll1.jaw_L, color='black', ls='--')
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