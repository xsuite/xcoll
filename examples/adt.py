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

# This example is a more realistic example of the ADT in the LHC. Each
# particle in a bunch gets the same kick (this is the default for a BlowUp).
# Note that this realistic behaviour shakes the beam around at every pass:
# the distribution becomes unmatched as it is no longer centered around
# the closed orbit. After the blow-up, the beam needs some relaxation time
# for the particles to decohere again.
# This is investigated in this example: the emittance is calculated by an
# emittance monitor (based on the covariance matrix of physical coordinates),
# and by a particles monitor that calculates the emittance based on normalised
# coordinates (which implicitly assumes a matched beam). In the resulting
# plots, the emittance based on normalised coordinates directly reflects its
# final value once the blow-up is stopped, while the emittance based on physical
# coordinates needs some relaxation time to reach this final value.

beam = 1
plane = 'V'
adt_turns = 1000
total_turns = 3500
num_particles = 5000
nemitt_x = 3.5e-6
nemitt_y = 2.5e-6


# Import a Run 3 LHC lattice without apertures
line = xt.Line.from_json(xc._pkg_root.parent / 'examples' / 'machines' / f'lhc_run3_b{beam}_no_aper.json')


# Create the ADT
pos = 'b5l4' if f'{beam}' == '1' and plane == 'H' else 'b5r4'
pos = 'b5l4' if f'{beam}' == '2' and plane == 'V' else pos
name = f'adtk{plane.lower()}.{pos}.b{beam}'
tank_start = f'adtk{plane.lower()}.{pos}.a.b{beam}'
tank_end   = f'adtk{plane.lower()}.{pos}.d.b{beam}'
adt_pos = 0.5*line.get_s_position(tank_start) + 0.5*line.get_s_position(tank_end)
adt = xc.BlowUp.install(line, name=name, at_s=adt_pos, need_apertures=False, plane=plane,
                        stop_at_turn=adt_turns)


# Add an emittance monitor
mon = xc.EmittanceMonitor.install(line, name="emittance monitor", at_s=adt_pos, stop_at_turn=total_turns)


# We also add a particles monitor, to calculate the emittance based on normalised coordinates.
# This will allow us to see how strongly the beam is mismatched.
mon2 = xt.ParticlesMonitor(start_at_turn=0, stop_at_turn=total_turns, num_particles=num_particles)
line.insert_element(element=mon2, name="particle monitor", at_s=adt_pos)


# Build the tracker and optimise
line.build_tracker()
line.optimize_for_tracking()
twiss = line.twiss()


# This will calibrate the ADT such that we push the full beam beyond 5 sigma over 1000 turns.
# In this example we won't blow up this aggressively (as the emittance calculation tends to
# get noisy at large emittances); hence we put the amplitude at 10%.
if plane == 'H':
    adt.calibrate_by_emittance(nemitt=nemitt_x, twiss=twiss)
else:
    adt.calibrate_by_emittance(nemitt=nemitt_y, twiss=twiss)
adt.amplitude = 0.1


# Generate a matched Gaussian bunch
part = xp.generate_matched_gaussian_bunch(num_particles=num_particles, total_intensity_particles=1.6e11,
                                          nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=7.55e-2, line=line)


# Move the tracker to a multi-core context
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))


# Activate the ADT
adt.activate()


# Track
line.track(part, num_turns=total_turns, with_progress=1)


# Get the normalised emittance from the particles monitor
part_norm = twiss.get_normalized_coordinates(mon2)
mon2_nemitt_x = np.array([np.mean(x[~np.isnan(x)])/2 for x in (part_norm.x_norm**2 + part_norm.px_norm**2).T])
mon2_nemitt_y = np.array([np.mean(y[~np.isnan(y)])/2 for y in (part_norm.y_norm**2 + part_norm.py_norm**2).T])
mon2_x_mean = np.mean(mon2.x, axis=0)
mon2_y_mean = np.mean(mon2.y, axis=0)
mon2_x_std =  np.std(mon2.x, axis=0)
mon2_y_std =  np.std(mon2.y, axis=0)
mon2_turns =  list(range(total_turns))


# Plot the results
_, ax = plt.subplots(figsize=(6,4))
ax.plot(mon.turns, 1.e6*mon.nemitt_x, label='H')
ax.plot(mon.turns, 1.e6*mon.nemitt_I, label='I')
ax.plot(mon2_turns, 1.e6 * mon2_nemitt_x * mon.beta0 * mon.gamma0, label='<xN^2 + pxN^2>/2')
ax.axvline(adt_turns, c='r', ls='--', label='stop blow-up')
ax.set_ylabel(r"$\epsilon\; [\mu\mathrm{m}]$")
ax.set_xlabel("Turn number")
ax.legend()
ax.set_title("Horizontal emittance growth by ADT blow-up in the LHC")
print(f"Total calculation time {time.time()-start_time}s")
plt.savefig("adt_horizontal_emittance.png", dpi=300)
plt.show()

_, ax = plt.subplots(figsize=(6,4))
ax.fill_between(mon2_turns, 1e3*mon2_x_mean + 1e3*mon2_x_std, 1e3*mon2_x_mean - 1e3*mon2_x_std, alpha=0.4)
ax.plot(mon2_turns, 1e3*mon2_x_mean, label=r'$<x>$')
ax.axvline(adt_turns, c='r', ls='--', label='stop blow-up')
ax.set_ylabel(r"normalised amplitude $[mm]$")
ax.set_xlabel("Turn number")
ax.legend()
ax.set_title("Average amplitude growth by ADT blow-up in the LHC")
plt.savefig("adt_horizontal_amplitude.png", dpi=300)
plt.show()

_, ax = plt.subplots(figsize=(6,4))
ax.plot(mon.turns, 1.e6*mon.nemitt_y, label='V')
ax.plot(mon.turns, 1.e6*mon.nemitt_II, label='II')
ax.plot(mon2_turns, 1.e6 * mon2_nemitt_y * mon.beta0 * mon.gamma0, label='<yN^2 + pyN^2>/2')
ax.axvline(adt_turns, c='r', ls='--', label='stop blow-up')
ax.set_ylabel(r"$\epsilon\; [\mu\mathrm{m}]$")
ax.set_xlabel("Turn number")
ax.legend()
ax.set_title("Vertical emittance growth by ADT blow-up in the LHC")
plt.savefig("adt_vertical_emittance.png", dpi=300)
plt.show()

_, ax = plt.subplots(figsize=(6,4))
ax.fill_between(mon2_turns, 1e3*mon2_y_mean + 1e3*mon2_y_std, 1e3*mon2_y_mean - 1e3*mon2_y_std, alpha=0.4)
ax.plot(mon2_turns, 100*mon2_y_mean, label=r'$<y>$')
ax.axvline(adt_turns, c='r', ls='--', label='stop blow-up')
ax.set_ylabel(r"normalised amplitude $[mm]$")
ax.set_xlabel("Turn number")
ax.legend()
ax.set_title("Average amplitude growth by ADT blow-up in the LHC")
plt.savefig("adt_vertical_amplitude.png", dpi=300)
plt.show()

_, ax = plt.subplots(figsize=(6,4))
ax.plot(mon.turns, 1.e2*mon.nemitt_zeta, label='L')
ax.plot(mon.turns, 1.e2*mon.nemitt_III, label='III')
ax.axvline(adt_turns, c='r', ls='--', label='stop blow-up')
ax.set_ylabel(r"$\epsilon\; [\mathrm{cm}]$")
ax.set_xlabel("Turn number")
ax.legend()
ax.set_title("Longitudinal emittance growth by ADT blow-up in the LHC")
plt.savefig("adt_longitudinal_emittance.png", dpi=300)
plt.show()
