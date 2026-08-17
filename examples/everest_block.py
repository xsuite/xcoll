# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import time
import numpy as np
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xcoll as xc

context = xo.ContextCpu()         # For CPU
# context = xo.ContextCupy()      # For CUDA GPUs
# context = xo.ContextPyopencl()  # For OpenCL GPUs

num_part = int(10e6)

block = xc.EverestBlock(length=1., material=xc.materials.Tungsten,
                        _context=context)

part = xt.Particles(x=np.zeros(num_part), energy0=450.e9, _context=context)

print(f"Tracking {num_part} particles through EverestBlock...   ", end='', flush=True)
t_start = time.time()
block.track(part)
print(f"Done in {time.time()-t_start:.2f} s")

x = context.nparray_from_context_array(part.x)
xp = context.nparray_from_context_array(part.px*part.rpp)

_, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].hist(x, bins=200, density=True)
ax[0].set_title('x distribution')
ax[1].hist(xp, bins=200, density=True)
ax[1].set_title('px distribution')
plt.tight_layout()
plt.show()
