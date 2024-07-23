# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc

beam = 1
plane = 'V'
num_turns = 100
num_particles = 5000
nemitt_x = 3.5e-6
nemitt_y = 2.5e-6


# Import a Run 3 LHC lattice without apertures
line = xt.Line.from_json(xc._pkg_root.parent / 'examples' / 'machines' / f'lhc_run3_b{beam}_no_aper.json')


# Create the ADT
adt = xc.BlowUp(plane=plane, amplitude=1)
pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
name = f'adtk{plane.lower()}.{pos}.b{beam}'
tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
adt.install(line, name=name, at_s=adt_pos, need_apertures=False)


# Add an emittance monitor
mon = xc.EmittanceMonitor(stop_at_turn=num_turns+1)
mon.set_beta_gamma_rel(line.particle_ref)
line.insert_element(element=mon, name="monitor", at_s=adt_pos)


# Build the tracker and optimise
line.build_tracker()
line.optimize_for_tracking()


# This will calibrate the ADT such that we gain ~ one emittance per turn.
# Note that this quickly explodes exponentially once the emittance becomes large.
if plane == 'H':
    adt.calibrate_by_emittance(nemitt=nemitt_x)
else:
    adt.calibrate_by_emittance(nemitt=nemitt_y)


# Generate a matched Gaussian bunch
part = xp.generate_matched_gaussian_bunch(num_particles=num_particles, total_intensity_particles=1.6e11,
                                          nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=7.55e-2, line=line)


# Move the tracker to a multi-core context
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))


# Activate the ADT
adt.activate()
adt.amplitude = 0.25


# Track
line.track(part, num_turns=num_turns, with_progress=1)


# Plot the result
_, ax = plt.subplots(figsize=(6,4))
t = list(range(num_turns+1))
ax.plot(t, 1.e6*mon.nemitt_x, label='H')
ax.plot(t, 1.e6*mon.nemitt_y, label='V')
ax.set_ylabel(r"$\epsilon\; [\mu\mathrm{m}]$")
ax.set_xlabel("turn")
ax.legend()
ax.set_title("Emittance growth by ADT blow-up in the LHC")
print(f"Total calculation time {time.time()-start_time}s")
plt.show()

_, ax = plt.subplots(figsize=(6,4))
ax.fill_between(t, mon.x_mean + mon.x_std, mon.x_mean - mon.x_std, alpha=0.2)
ax.plot(t, mon.x_mean, label=r'$<x_N>$')
ax.fill_between(t, mon.y_mean + mon.y_std, mon.y_mean - mon.y_std, alpha=0.2)
ax.plot(t, mon.y_mean, label=r'$<y_N>$')
ax.set_ylabel(r"normalised amplitude $[\sigma]$")
ax.set_xlabel("turn")
ax.legend()
ax.set_title("Average amplitude growth by ADT blow-up in the LHC")
plt.show()
