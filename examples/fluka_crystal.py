# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


# Make a context and get a buffer
context = xo.ContextCpu()         # For CPU
# context = xo.ContextCupy()      # For CUDA GPUs
# context = xo.ContextPyopencl()  # For OpenCL GPUs

num_part = int(10_000)

coll = xc.FlukaCrystal(length=0.002, material='si', bending_angle=149e-6,
                       width=0.002, height=0.05, side='+', jaw=1e-3,  # Hack to make it work for now; jaw should be almost zero
                         _context=context)

# Connect to FLUKA
xc.fluka.engine.particle_ref = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=4e11)
xc.fluka.engine.capacity = 2*num_part
xc.fluka.engine.seed = 5656565
xc.fluka.engine.start(elements=coll, clean=True, verbose=False)

x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init, particle_ref=xc.fluka.engine.particle_ref, _context=context, _capacity=xc.fluka.engine.capacity)
part_init = part.copy()

coll.track(part)

xc.fluka.engine.stop()

# Sort particles to be able to compare to part_init
part.sort(interleave_lost_particles=True)

# Select only surviving particles and only within the window of interest
mask = (part.state > 0 ) & ( part.px - part_init.px < 250.e-6) & ( part.px - part_init.px > -50.e-6)

plt.figure(figsize=(15,10))
plt.hist2d(part_init.px[mask]*1.e6, part.px[mask]*1.e6 - part_init.px[mask]*1.e6, 500, norm=mpl.colors.LogNorm())
plt.xlim(-30, 180)
plt.ylim(-55, 205)
plt.ylabel(r'$\Delta\theta$ [$\mu$rad]')
plt.xlabel(r'$\theta_{in}$ [$\mu$rad]')
plt.tight_layout()
plt.show()

