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


# Make a context and get a buffer
context = xo.ContextCpu(omp_num_threads='auto') # For CPU
# context = xo.ContextCupy()                    # For CUDA GPUs
# context = xo.ContextPyopencl()                # For OpenCL GPUs

num_part = int(5e6)
particle_ref = xt.Particles('proton', p0c=4e11)

radii = [100, 50, 10, 5, 3, 2, 1.5, 1, 0.5]

for R in radii:
    coll = xc.EverestPreciseCrystal(length=0.002, material=xc.materials.Silicon, bending_radius=R,
                                    width=0.002, height=0.05, side='+', miscut=0., lattice='strip', jaw=0.001,
                                    _context=context)
    ang = coll.bending_angle * 1e6
    print(f"R={R} m, bending angle = {ang:.2f} urad")

    x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
    px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                            particle_ref=particle_ref, _context=context)
    part_init = part.copy()

    print(f"Tracking {num_part} particles...     ", end='', flush=True)
    start = time.time()
    coll.track(part)
    print(f"Done in {time.time() - start:.2f} seconds.")

    # Sort particles to be able to compare to part_init
    part.sort(interleave_lost_particles=True)

    # Select only surviving particles and only within the window of interest
    mask = (part.state > 0 ) & ( part.px - part_init.px < (ang+55)*1.e-6) & ( part.px - part_init.px > -50.e-6)

    fig, ax = plt.subplots(1, 2, figsize=(15,7))
    ax[0].hist2d(part_init.px[mask]*1.e6, part.px[mask]*1.e6 - part_init.px[mask]*1.e6, 500, norm=mpl.colors.LogNorm())
    ax[0].set_xlim(-30, 180)
    ax[0].set_ylim(-55, ang+55)
    ax[0].set_ylabel(r'$\Delta\theta$ [$\mu$rad]')
    ax[0].set_xlabel(r'$\theta_{in}$ [$\mu$rad]')

    # Zoom on channelling region

    x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
    px_init  = np.random.uniform(low=-13.e-6, high=13.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                            particle_ref=particle_ref, _context=context)
    part_init = part.copy()

    print(f"Tracking {num_part} particles (all in channel)...     ", end='', flush=True)
    start = time.time()
    coll.track(part)
    print(f"Done in {time.time() - start:.2f} seconds.")

    # Sort particles to be able to compare to part_init
    part.sort(interleave_lost_particles=True)

    # Select only surviving particles and only within the window of interest
    mask = (part.state > 0 ) & ( part.px < (ang+13)*1.e-6) & ( part.px > (ang-13)*1.e-6)

    ax[1].hist2d(part_init.px[mask]*1.e6, part.px[mask]*1.e6 , 500, norm=mpl.colors.LogNorm())
    ax[1].set_xlim(-13, 13)
    ax[1].set_ylim(ang-13, ang+13)
    ax[1].set_ylabel(r'$\theta_{out}$ [$\mu$rad]')
    ax[1].set_xlabel(r'$\theta_{in}$ [$\mu$rad]')

    plt.tight_layout()
    plt.savefig(f'plots/everest_precise_crystal_{R}.png', dpi=300)

plt.show()
