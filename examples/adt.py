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
n_turns = 100
num_particles = 5000
nemitt_x = 3.5e-6
nemitt_y = 3.5e-6


# Helper function to calculate the emittance
def calculate_nemitt(part):
    cov_x = np.cov(part.x, part.px)
    cov_y = np.cov(part.y, part.py)
    nemitt_x = part.beta0[0]*part.gamma0[0]*np.sqrt(cov_x[0,0]*cov_x[1,1]-cov_x[1,0]*cov_x[0,1])
    nemitt_y = part.beta0[0]*part.gamma0[0]*np.sqrt(cov_y[0,0]*cov_y[1,1]-cov_y[1,0]*cov_y[0,1])
    return nemitt_x, nemitt_y


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


# Log the initial emittance and normalised amplitudes
ex, ey = calculate_nemitt(part)
ex = [ex]
ey = [ey]
tw = line.twiss()
part_norm = tw.get_normalized_coordinates(part, nemitt_x=nemitt_x, nemitt_y=nemitt_y)
x_norm = np.sqrt(part_norm.x_norm**2 + part_norm.px_norm**2)
x_norm_mean = [np.mean(x_norm)]
x_norm_std  = [np.std(x_norm)]
y_norm = np.sqrt(part_norm.y_norm**2 + part_norm.py_norm**2)
y_norm_mean = [np.mean(y_norm)]
y_norm_std  = [np.std(y_norm)]


# Move the tracker to a multi-core context
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))


# Activate the ADT
adt.activate()
adt.amplitude = 0.25


# Track and store emittance and normalised amplitude at every turn
for _ in range(n_turns):
    line.track(part)
    this_ex, this_ey = calculate_nemitt(part)
    ex.append(this_ex)
    ey.append(this_ey)
    part_norm = tw.get_normalized_coordinates(part, nemitt_x=nemitt_x, nemitt_y=nemitt_y)
    x_norm = np.sqrt(part_norm.x_norm**2 + part_norm.px_norm**2)
    x_norm_mean.append(np.mean(x_norm))
    x_norm_std.append(np.std(x_norm))
    y_norm = np.sqrt(part_norm.y_norm**2 + part_norm.py_norm**2)
    y_norm_mean.append(np.mean(y_norm))
    y_norm_std.append(np.std(y_norm))


print(f"Total calculation time {time.time()-start_time}s")


# Plot the result
_, ax = plt.subplots(figsize=(6,4))
s = list(range(n_turns+1))
ax.plot(s, 1.e6*np.array(ex), label='H')
ax.plot(s, 1.e6*np.array(ey), label='V')
ax.set_ylabel(r"$\epsilon\; [\mu\mathrm{m}]$")
ax.set_xlabel("turn")
ax.legend()
ax.set_title("Emittance growth by ADT blow-up in the LHC")
plt.show()

_, ax = plt.subplots(figsize=(6,4))
ax.fill_between(s, np.array(x_norm_mean) + np.array(x_norm_std), np.array(x_norm_mean) - np.array(x_norm_std), alpha=0.2)
ax.plot(s, np.array(x_norm_mean), label=r'$<x_N>$')
ax.fill_between(s, np.array(y_norm_mean) + np.array(y_norm_std), np.array(y_norm_mean) - np.array(y_norm_std), alpha=0.2)
ax.plot(s, np.array(y_norm_mean), label=r'$<y_N>$')
ax.set_ylabel(r"normalised amplitude $[\sigma]$")
ax.set_xlabel("turn")
ax.legend()
ax.set_title("Average amplitude growth by ADT blow-up in the LHC")
plt.show()
