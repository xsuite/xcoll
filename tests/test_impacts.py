# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt

# import xobjects as xo
# import xtrack as xt
# import xpart as xp
# import xcoll as xc


# num_part     = int(1e6)
# impacts_size = int(2e7) # 4-5 GB of memory


# # Make a context
# context = xo.ContextCpu()         # For CPU
# # context = xo.ContextCupy()      # For CUDA GPUs
# # context = xo.ContextPyopencl()  # For OpenCL GPUs


# # Make a block of material
# coll = xc.EverestBlock(length=0.2, material=xc.materials.Tungsten, _context=context)


# # Create initial particles (with zero angles)
# x_init    = np.random.normal(loc=0., scale=1e-3, size=num_part)
# y_init    = np.random.normal(loc=0., scale=1e-3, size=num_part)
# particles = xp.Particles(x=x_init, y=y_init, p0c=450e9, _context=context)


# # Initialise impacts logging
# io_buffer = xt.new_io_buffer(capacity=impacts_size)
# impacts = xt.start_internal_logging(elements=[coll], io_buffer=io_buffer, capacity=io_buffer.capacity)


# # Track !
# coll.track(particles)

# raise ValueError("Example not finished")
