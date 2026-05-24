# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import gc
import time
import numpy as np
from pathlib import Path
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

total_particles = int(250e6)
batch_size      = int(5e6)
particle_ref = xt.Particles('proton', p0c=4e11)

coll = xc.EverestPreciseCrystal(length=0.002, material=xc.materials.Silicon, bending_angle=149e-6,
                                width=0.002, height=0.05, side='+', miscut=0., lattice='strip', jaw=0.001,
                                _context=context)

def to_np(context, arr):
    # Portable CPU/GPU conversion
    return context.nparray_from_context_array(arr)


def accumulate_hist(
    *,
    total_particles,
    batch_size,
    make_initial_conditions,
    xedges,
    yedges,
    y_quantity,
    mask_function,
    label,
):
    H = np.zeros((len(xedges) - 1, len(yedges) - 1), dtype=np.float64)

    done = 0
    t_tot = 0
    while done < total_particles:
        n = min(batch_size, total_particles - done)

        x_init, px_init, y_init, py_init = make_initial_conditions(n)

        part = xp.build_particles(
            x=x_init, px=px_init, y=y_init, py=py_init,
            particle_ref=particle_ref,
            _context=context,
        )

        # Free arrays that are already copied into `part`
        del x_init, y_init, py_init

        print(f"{label}: tracking batch {int((done+n)/1e6)}M/{int(total_particles/1e6)}M ... ", end="", flush=True)
        start = time.time()
        coll.track(part)
        t_track = time.time() - start
        print(f"{t_track:.2f} s")
        t_tot += t_track

        part.sort(interleave_lost_particles=True)

        ids = to_np(context, part.particle_id).astype(np.int64)
        px0 = px_init[ids]
        px1 = to_np(context, part.px)
        state = to_np(context, part.state)

        mask = mask_function(px0, px1, state)
        y = y_quantity(px0, px1)

        H += np.histogram2d(
            px0[mask] * 1e6,
            y[mask] * 1e6,
            bins=(xedges, yedges),
        )[0]

        done += n

        del part, px_init, px0, px1, state, ids, mask, y
        gc.collect()

    return H, t_tot



# --------------------------------------------------------------------
# Full angular scan
# --------------------------------------------------------------------

xedges0 = np.linspace(-30, 180, 501)
yedges0 = np.linspace(-55, 205, 501)

def make_full(n):
    return (
        np.random.normal(loc=1.5e-3, scale=75.e-6, size=n),
        np.random.uniform(low=-50.e-6, high=250.e-6, size=n),
        np.random.normal(loc=0., scale=1e-3, size=n),
        np.random.normal(loc=0., scale=5.e-6, size=n),
    )

H0, t_tot0 = accumulate_hist(
    total_particles=total_particles,
    batch_size=batch_size,
    make_initial_conditions=make_full,
    xedges=xedges0,
    yedges=yedges0,
    y_quantity=lambda px0, px1: px1 - px0,
    mask_function=lambda px0, px1, state:
        (state > 0) & (px1 - px0 < 250.e-6) & (px1 - px0 > -50.e-6),
    label="Full scan",
)

# --------------------------------------------------------------------
# Channelling zoom
# --------------------------------------------------------------------

xedges1 = np.linspace(-13, 13, 501)
yedges1 = np.linspace(136, 162, 501)

def make_channel(n):
    return (
        np.random.normal(loc=1.5e-3, scale=75.e-6, size=n),
        np.random.uniform(low=-13.e-6, high=13.e-6, size=n),
        np.random.normal(loc=0., scale=1e-3, size=n),
        np.random.normal(loc=0., scale=5.e-6, size=n),
    )

H1, t_tot1 = accumulate_hist(
    total_particles=total_particles,
    batch_size=batch_size,
    make_initial_conditions=make_channel,
    xedges=xedges1,
    yedges=yedges1,
    y_quantity=lambda px0, px1: px1,
    mask_function=lambda px0, px1, state:
        (state > 0) & (px1 < 162.e-6) & (px1 > 136.e-6),
    label="Channel zoom",
)

print()
print(f"Total time for full scan: {t_tot0:.2f} s")
print(f"Total time for channel zoom: {t_tot1:.2f} s")


fig, ax = plt.subplots(1, 2, figsize=(15, 7))

pcm0 = ax[0].pcolormesh(
    xedges0, yedges0,
    np.ma.masked_where(H0.T == 0, H0.T),
    norm=mpl.colors.LogNorm(),
)
ax[0].set_xlim(-30, 180)
ax[0].set_ylim(-55, 205)
ax[0].set_ylabel(r'$\Delta\theta$ [$\mu$rad]')
ax[0].set_xlabel(r'$\theta_{in}$ [$\mu$rad]')

pcm1 = ax[1].pcolormesh(
    xedges1, yedges1,
    np.ma.masked_where(H1.T == 0, H1.T),
    norm=mpl.colors.LogNorm(),
)
ax[1].set_xlim(-13, 13)
ax[1].set_ylim(136, 162)
ax[1].set_ylabel(r'$\theta_{out}$ [$\mu$rad]')
ax[1].set_xlabel(r'$\theta_{in}$ [$\mu$rad]')

plt.tight_layout()
Path("plots").mkdir(exist_ok=True)
plt.savefig("plots/everest_precise_crystal_highres.png", dpi=300)
plt.show()
