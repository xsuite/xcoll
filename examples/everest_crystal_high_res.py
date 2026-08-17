# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


context = xo.ContextCpu()         # For CPU
# context = xo.ContextCupy()      # For CUDA GPUs
# context = xo.ContextPyopencl()  # For OpenCL GPUs

num_part = int(1e9)
chunk_size = int(10e6)
px_in_lim = (-50e-6, 250e-6)
px_out_lim = (-50e-6, 250e-6)
print_n_chunks = 10
nbins = 500

particle_ref = xt.Particles('proton', p0c=4e11)
coll = xc.EverestCrystal(length=0.002, material=xc.materials.Silicon, bending_angle=149e-6,
                         width=0.002, height=0.05, side='+', miscut=0., lattice='strip', jaw=0.001,
                         _context=context)

# Tracking init
xnp = context.nplike_lib
steps = (num_part - 1) // chunk_size + 1
total_surv = 0

# Histogram init
xrange = (px_in_lim[0] * 1.e6, px_in_lim[1] * 1.e6)
yrange = (px_out_lim[0] * 1.e6, px_out_lim[1] * 1.e6)
xedges = np.linspace(*xrange, nbins + 1)
yedges = np.linspace(*yrange, nbins + 1)
H = xnp.zeros((nbins, nbins), dtype=xnp.uint64)

for i in range(steps):
    if i == steps - 1:
        # Last step
        chunk_size = num_part - i*chunk_size
    if i % print_n_chunks == 0:
        ch = f"{i+1}-{i+print_n_chunks}"
        print(f"Tracking chunks {ch:>7}/{steps} ("
              f"{print_n_chunks*chunk_size} particles)...   ",
              end='', flush=True)
        t_part = 0
        t_track = 0
        t_sort = 0
        t_post = 0

    t_start = time.time()
    x_init  = np.random.normal(loc=1.5e-3, scale=75.e-6, size=chunk_size)
    px_init = np.random.uniform(low=px_in_lim[0], high=px_in_lim[1], size=chunk_size)
    y_init  = np.random.normal(loc=0., scale=1e-3, size=chunk_size)
    py_init = np.random.normal(loc=0., scale=5.e-6, size=chunk_size)
    part = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                              particle_ref=particle_ref, _context=context)
    px_init = part.px.copy()
    t_part += time.time() - t_start

    t_start = time.time()
    coll.track(part)
    t_track += time.time() - t_start

    # Sort particles to be able to compare to part_init
    # Do this on the original context to save computation time
    t_start = time.time()
    idx = xnp.argsort(part.particle_id)
    px = part.px[idx]
    state = part.state[idx]
    t_sort += time.time() - t_start

    # Select only surviving particles and only within the window of interest
    t_start = time.time()
    mask = (state > 0 ) & (px - px_init > px_out_lim[0]) & (px - px_init < px_out_lim[1])
    num_surv = mask.sum().tolist()
    total_surv += num_surv
    H_chunk, _, _ = xnp.histogram2d(
                        px_init[mask]*1.e6, (px[mask]-px_init[mask]) * 1.e6,
                        bins=nbins, range=(xrange, yrange),
                    )
    H += H_chunk.astype(xnp.uint64)
    del px, px_init, state, mask, part, idx, H_chunk # To force-free memory
    t_post += time.time() - t_start

    if i % print_n_chunks == print_n_chunks - 1:
        print(f"Done in {t_part:.2f}s (creating particles) + "
              f"{t_track:.2f}s (tracking) + {t_sort:.2f}s (sorting) + "
              f"{t_post:.2f}s (post-processing).")

print(f"Out of {num_part} particles, {total_surv} survived.")
H_cpu = context.nparray_from_context_array(H)
plt.figure(figsize=(12, 8))
plt.pcolormesh(xedges, yedges, H_cpu.T, norm=mpl.colors.LogNorm())
plt.xlim(*xrange)
plt.ylim(*yrange)
plt.ylabel(r'$\Delta\theta$ [$\mu$rad]')
plt.xlabel(r'$\theta_{in}$ [$\mu$rad]')
plt.tight_layout()
plt.savefig('everest_crystal_high_res.png', dpi=300)
plt.show()
