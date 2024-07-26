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
adt_turns = 1000
total_turns = 3500
num_particles = 5000
nemitt_x = 3.5e-6
nemitt_y = 3.5e-6


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
mon = xc.EmittanceMonitor(stop_at_turn=total_turns)
mon.set_beta0_gamma0(line.particle_ref)
line.insert_element(element=mon, name="emittance monitor", at_s=adt_pos)


# Add a particles monitor
mon2 = xt.ParticlesMonitor(start_at_turn=0, stop_at_turn=total_turns, num_particles=num_particles)
line.insert_element(element=mon2, name="particle monitor", at_s=adt_pos)


# Install BlackAbsorbers
xc.install_elements(line, names=['tcp.c6l7.b1', 'tcp.d6l7.b1', 'tcp.b6l7.b1'], elements=[
        xc.BlackAbsorber(length=1, gap=5, angle=0),
        xc.BlackAbsorber(length=1, gap=5, angle=90),
        xc.BlackAbsorber(length=1, gap=5, angle=135)])


# Build the tracker and optimise
line.build_tracker()
line.optimize_for_tracking()


# This will calibrate the ADT such that we gain ~ one emittance per turn.
# Note that this quickly explodes exponentially once the emittance becomes large.
if plane == 'H':
    adt.calibrate_by_emittance(nemitt=nemitt_x)
else:
    adt.calibrate_by_emittance(nemitt=nemitt_y)
twiss = line.twiss()
xc.assign_optics_to_collimators(line, nemitt_x=3.5e-6, nemitt_y=3.5e-6, twiss=twiss)


# Generate a matched Gaussian bunch
part = xp.generate_matched_gaussian_bunch(num_particles=num_particles, total_intensity_particles=1.6e11,
                                          nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=7.55e-2, line=line)
part._init_random_number_generator(seeds=23*np.ones(num_particles, dtype=np.int64))


# Move the tracker to a multi-core context
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))


# Activate the ADT
adt.activate()
adt.amplitude = 0.025
print(adt._kick_rms)
xc.enable_scattering(line)


# Track
line.track(part, num_turns=adt_turns, with_progress=1)
at_ele, counts = np.unique(part.at_element, return_counts=True)
for el,c in zip(at_ele, counts):
    print(f"{line.element_names[el]} {c}")


# Plot the result
_, ax = plt.subplots(figsize=(6,4))
t = list(range(adt_turns))
ax.plot(t, 1.e6*mon.nemitt_x, label='H')
ax.plot(t, 1.e6*mon.nemitt_y, label='V')
ax.plot(t, 1.e6*mon.nemitt_I, label='I')
ax.plot(t, 1.e6*mon.nemitt_II, label='II')
ax.set_ylabel(r"$\epsilon\; [\mu\mathrm{m}]$")
ax.set_xlabel("turn")
ax.legend()
ax.set_title("Emittance growth by ADT blow-up in the LHC")
print(f"Total calculation time {time.time()-start_time}s")
plt.show()


# Track
adt.deactivate()
line.track(part, num_turns=total_turns-adt_turns, with_progress=1)
at_ele, counts = np.unique(part.at_element, return_counts=True)
for el,c in zip(at_ele, counts):
    print(f"{line.element_names[el]} {c}")
part_norm = twiss.get_normalized_coordinates(mon2, nemitt_x=3.5e-6, nemitt_y=3.6e-6)
part_norm.to_pandas().to_csv("adt_norm_particles.csv", index=False)


# Plot the result
_, ax = plt.subplots(figsize=(6,4))
t = list(range(total_turns))
ax.plot(t, 1.e6*mon.nemitt_x, label='H')
ax.plot(t, 1.e6*mon.nemitt_y, label='V')
ax.plot(t, 1.e6*mon.nemitt_I, label='I')
ax.plot(t, 1.e6*mon.nemitt_II, label='II')
ax.axvline(adt_turns, c='r', ls='--', label='stop blow-up')
ax.set_ylabel(r"$\epsilon\; [\mu\mathrm{m}]$")
ax.set_xlabel("turn")
ax.legend()
ax.set_title("Emittance growth by ADT blow-up in the LHC")
print(f"Total calculation time {time.time()-start_time}s")
plt.savefig("adt_transverse_emittances.png", dpi=300)
plt.show()

_, ax = plt.subplots(figsize=(6,4))
ax.fill_between(t, 1e3*mon.x_mean + 1e3*mon.x_std, 1e3*mon.x_mean - 1e3*mon.x_std, alpha=0.4)
ax.plot(t, 1e3*mon.x_mean, label=r'$<x>$')
ax.fill_between(t, 1e3*mon.y_mean + 1e3*mon.y_std, 1e3*mon.y_mean - 1e3*mon.y_std, alpha=0.4)
ax.plot(t, 100*mon.y_mean, label=r'$<y>$')
ax.axvline(adt_turns, c='r', ls='--', label='stop blow-up')
ax.set_ylabel(r"normalised amplitude $[mm]$")
ax.set_xlabel("turn")
ax.legend()
ax.set_title("Average amplitude growth by ADT blow-up in the LHC")
plt.savefig("adt_transverse_amplitudes.png", dpi=300)
plt.show()

_, ax = plt.subplots(figsize=(6,4))
t = list(range(total_turns))
ax.plot(t, 1.e2*mon.nemitt_zeta, label='L')
ax.plot(t, 1.e2*mon.nemitt_III, label='III')
ax.axvline(adt_turns, c='r', ls='--', label='stop blow-up')
ax.set_ylabel(r"$\epsilon\; [\mathrm{cm}]$")
ax.set_xlabel("turn")
ax.legend()
ax.set_title("Longitudinal emittance growth by ADT blow-up in the LHC")
plt.savefig("adt_longitudinal_emittances.png", dpi=300)
plt.show()
