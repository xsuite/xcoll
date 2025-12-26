# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


num_part = int(2_000)
capacity = 2*num_part
particle_ref = xt.Particles('proton', p0c=4e11)

coll = xc.FlukaCrystal(length=0.002, material=xc.materials.SiliconCrystal, bending_angle=149e-6,
                       width=0.002, height=0.05, side='+', jaw=0.001)

# Connect to FLUKA
xc.fluka.engine.particle_ref = particle_ref
xc.fluka.engine.capacity = capacity
xc.fluka.engine.seed = 5656565
xc.fluka.engine.start(elements=coll, clean=False, verbose=False)

x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                          particle_ref=xc.fluka.engine.particle_ref,
                          _capacity=xc.fluka.engine.capacity)
part_init = part.copy()

print(f"Tracking {num_part} particles (FLUKA)...     ", end='', flush=True)
start = time.time()
coll.track(part)
print(f"Done in {time.time() - start:.2f} seconds.")

# Stop the FLUKA server
xc.fluka.engine.stop(clean=True)

# Sort particles to be able to compare to part_init
part.sort(interleave_lost_particles=True)

# Select only surviving particles and only within the window of interest
mask = (part.state > 0 ) & ( part.px - part_init.px < 250.e-6) & ( part.px - part_init.px > -50.e-6)

plt.figure(figsize=(12,8))
plt.hist2d(part_init.px[mask]*1.e6, part.px[mask]*1.e6 - part_init.px[mask]*1.e6, 500, norm=mpl.colors.LogNorm())
plt.xlim(-30, 180)
plt.ylim(-55, 205)
plt.ylabel(r'$\Delta\theta$ [$\mu$rad]')
plt.xlabel(r'$\theta_{in}$ [$\mu$rad]')
plt.tight_layout()
plt.show()

