import time
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc


num_part = 100000
_capacity = 2*num_part
particle_ref = xp.Particles('proton', p0c=7e12)

everest_coll = xc.EverestCollimator(length=1.0, material=xc.materials.MolybdenumGraphite, jaw=[0.01, -0.01])
# fluka_coll = xc.FlukaCollimator(length=1.0, jaw=[0.01, -0.01])
g4_coll = xc.Geant4Collimator(length=1.0, material='mogr', jaw=[0.1, -0.1], geant4_id='g4_coll')


# Plotting function
# -----------------

dri = xt.Drift(length=everest_coll.length)
def plot_scatters(part_init, part, outfile, show_colorbar=True):
    devmin = 1e-6
    devmax = 1
    part_dri = part_init.copy()
    dri.track(part_dri)
    part_dri.sort(interleave_lost_particles=True)
    part.sort(interleave_lost_particles=True)
    mask = part.state > 0
    dev = 1.e3*np.sqrt((part.x[mask] - part_dri.x[part.parent_particle_id[mask]])**2
                     + (part.y[mask] - part_dri.y[part.parent_particle_id[mask]])**2)
    dev[dev<devmin] = devmin
    if len(dev[dev>devmax]) > 0:
        print(f'Warning: {len(dev[dev>devmax])} particles have deviation > {devmax} mm')

    if show_colorbar:
        fig = plt.figure(figsize=(7.11, 3.5))
    else:
        fig = plt.figure(figsize=(6, 3.5))
    plt.scatter(1.e3*part_init.x[part.parent_particle_id[mask]],
                1.e3*part_init.y[part.parent_particle_id[mask]],
                c=dev, s=0.25, cmap='viridis',norm=mcolors.LogNorm(vmin=devmin, vmax=devmax))
    plt.vlines(10, -5, 5, colors='r')
    plt.xlim(8.5, 11.5)
    plt.ylim(-5, 5)
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    if show_colorbar:
        plt.colorbar(label='Deviation from drift [mm]')
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)


# Generate particles
# ------------------

y = np.random.normal(loc=0, scale=1e-3, size=num_part)
py = np.random.normal(loc=0, scale=1e-3, size=num_part)
delta = np.random.normal(loc=0, scale=3e-4, size=num_part)

part_init = []
x0 = 9.25e-3
x = np.random.normal(loc=x0, scale=1e-4, size=num_part)
part_init.append(xp.build_particles(particle_ref=particle_ref, x=x, px=0, y=y, py=py, delta=delta, _capacity=_capacity))

x0 = 0.01
x = np.random.normal(loc=x0, scale=1e-4, size=num_part)
part_init.append(xp.build_particles(particle_ref=particle_ref, x=x, px=0, y=y, py=py, delta=delta, _capacity=_capacity))

x0 = 0.01075
x = np.random.normal(loc=x0, scale=1e-4, size=num_part)
part_init.append(xp.build_particles(particle_ref=particle_ref, x=x, px=0, y=y, py=py, delta=delta, _capacity=_capacity))


# # Fluka
# # -----

# xc.FlukaEngine.start(elements=fluka_coll, names='tcp.c6l7.b1', debug_level=1, _capacity=_capacity, cwd='temp_scratch')
# xc.FlukaEngine.set_particle_ref(particle_ref=particle_ref)
# part_fluka = []
# for this_part_init in part_init:
#     this_part = this_part_init.copy()
#     t_start = time.time()
#     fluka_coll.track(this_part)
#     t_end = time.time()

#     print(f"FLUKA execution time: {t_end - t_start} s")
#     print(f"Survived in FLUKA: {len(np.unique(this_part.parent_particle_id[this_part.state>0]))}/{num_part}")
#     this_part.sort(interleave_lost_particles=True)
#     part_fluka.append(this_part)

# xc.FlukaEngine.stop()

# for idx in range(3):
#     plot_scatters(part_init[idx], part_fluka[idx], f'fluka_collimator_{idx}.png', show_colorbar=True)


# Everest
# -------

# Pre-run for compilation
this_part = part_init[0].copy()
everest_coll.track(this_part)

part_everest = []
for this_part_init in part_init:
    this_part = this_part_init.copy()
    t_start = time.time()
    everest_coll.track(this_part)
    t_end = time.time()

    print(f"Everest execution time: {t_end - t_start} s")
    print(f"Survived in Everest: {len(np.unique(this_part.parent_particle_id[this_part.state>0]))}/{num_part}")
    this_part.sort(interleave_lost_particles=True)
    part_everest.append(this_part)

for idx in range(3):
    plot_scatters(part_init[idx], part_everest[idx], f'everest_collimator_{idx}.png', show_colorbar=True)


# Geant4
# ------

xc.Geant4Engine.start(elements=g4_coll, seed=1993, bdsim_config_file='settings_protons.gmad', particle_ref=particle_ref)

part_g4 = []
for this_part_init in part_init:
    this_part = this_part_init.copy()
    t_start = time.time()
    g4_coll.track(this_part)
    t_end = time.time()

    print(f"Geant4 execution time: {t_end - t_start} s")
    print(f"Survived in Geant4: {len(np.unique(this_part.parent_particle_id[this_part.state>0]))}/{num_part}")
    this_part.sort(interleave_lost_particles=True)
    part_g4.append(this_part)

xc.Geant4Engine.stop()

for idx in range(3):
    plot_scatters(part_init[idx], part_g4[idx], f'g4_collimator_{idx}.png', show_colorbar=True)
