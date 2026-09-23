# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

import xobjects as xo
import xtrack as xt
import xcoll as xc

######################################################
# Constants
######################################################
KB = sc.Boltzmann
T_ROOM = 293.15  # K

######################################################
# Beam parameters
######################################################
nemitt_x = 1e-5
nemitt_y = 1e-7

sigma_z = 4e-3
sigma_delta = 1e-3

bunch_intensity = 4e9

######################################################
# Build a toy ring
######################################################
lbend = 3
angle = np.pi / 2

lquad = 0.3
k1qf = 0.1
k1qd = 0.7

env = xt.Environment()

line = env.new_line(components=[
    env.new('mqf.1', xt.Quadrupole, length=lquad, k1=k1qf),
    env.new('d1.1',  xt.Drift, length=1),
    env.new('mb1.1', xt.Bend, length=lbend, angle=angle),
    env.new('d2.1',  xt.Drift, length=1),

    env.new('mqd.1', xt.Quadrupole, length=lquad, k1=-k1qd),
    env.new('d3.1',  xt.Drift, length=1),
    env.new('mb2.1', xt.Bend, length=lbend, angle=angle),
    env.new('d4.1',  xt.Drift, length=1),

    env.new('mqf.2', xt.Quadrupole, length=lquad, k1=k1qf),
    env.new('d1.2',  xt.Drift, length=1),
    env.new('mb1.2', xt.Bend, length=lbend, angle=angle),
    env.new('d2.2',  xt.Drift, length=1),

    env.new('mqd.2', xt.Quadrupole, length=lquad, k1=-k1qd),
    env.new('d3.2',  xt.Drift, length=1),
    env.new('mb2.2', xt.Bend, length=lbend, angle=angle),
    env.new('d4.2',  xt.Drift, length=1),
])

line.set_particle_ref('electron', p0c=1e9)
line.configure_bend_model(core='full', edge=None)

######################################################
# Insert beam-gas scattering centers
######################################################
# Each BeamGasScattering element represents the lattice section between the
# previous scattering center (or the start of the line) and itself, so the
# last one must sit at the end of the line to cover the whole circumference.
circumference = line.get_length()
n_beamgas = 16

placements = []
for ii, ss in enumerate(np.linspace(0, circumference, n_beamgas + 1)[1:]):
    beamgas_name = f'BeamGasScattering.{ii}'
    env.elements[beamgas_name] = xc.BeamGasScattering()
    placements.append(env.place(beamgas_name, at=ss))

line.insert(placements)

######################################################
# Install apertures
######################################################
tab = line.get_table()
needs_aperture = tab.rows.match_not(
    element_type='Drift.*|Marker|').name

aper_size = 0.040  # m

env.new('aper', xt.LimitRect,
        min_x=-aper_size, max_x=aper_size,
        min_y=-aper_size, max_y=aper_size)

placements = []
for nn in needs_aperture:
    env.new(f'{nn}_aper_entry', 'aper')
    env.new(f'{nn}_aper_exit', 'aper')
    placements.append(env.place(f'{nn}_aper_entry', at=f'{nn}@start'))
    placements.append(env.place(f'{nn}_aper_exit', at=f'{nn}@end'))

line.insert(placements)

######################################################
# Define an example flat pressure profile
######################################################
tab = line.get_table()
tt_beamgas = tab.rows[tab.element_type == 'BeamGasScattering']

# N2 at 1e-7 mbar and room temperature. The gas density table holds the
# *atomic* density of each species, hence the factor 2 for the N2 molecule.
_mbar_to_pascal = 1e2
pressure_mbar = 1e-7
pressure_pascal = pressure_mbar * _mbar_to_pascal

atomic_density = 2 * pressure_pascal / (KB * T_ROOM)

gas_density = xt.Table({
    'name': tt_beamgas.name,
    's': tt_beamgas.s,
    'N': np.ones(len(tt_beamgas.name)) * atomic_density,
})

######################################################
# Beam-gas simulation
######################################################
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads='auto'))

# The facade builds the study and initialises it in one call. Equivalently,
# construct xc.BeamGasStudy(line=line, ...) and call initialise_beamgas().
beamgas = line.xcoll.beamgas_configure(
    gas_density=gas_density,
    process='coulomb',
    # Only the large angles can drive a particle into the aperture; restricting
    # the generated range is the main variance-reduction knob of the study.
    coulomb_theta=(8e-3, 20e-3),
    # process='brems',
    # brems_energy_cut=1e6,
    nemitt_x=nemitt_x,
    nemitt_y=nemitt_y,
    sigma_z=sigma_z,
    sigma_delta=sigma_delta,
    bunch_intensity=bunch_intensity,
    n_scattering_events=1000,
    seed=1997,
    method='4d',
)

print(beamgas.local_rates())

result = beamgas.run(
    track=True,
    n_turns=100,
    keep_particles=True,
    with_progress=1,
)

print(f'Beam-gas interaction rate: {result.rate_scattering*1e-3:.3f} kHz')
print(f'Beam-gas loss rate:        {result.rate_tracking*1e-3:.3f} kHz')
print(f'Beam-gas lifetime:         {result.lifetime_tracking/60:.2f} min')

######################################################
# Optional: refine loss locations
######################################################
loss_loc_refinement = xt.LossLocationRefinement(
    line,
    n_theta=360,  # Angular resolution of the polygonal aperture approximation
    r_max=0.5,    # Maximum transverse aperture in m
    dr=50e-6,     # Transverse loss refinement accuracy [m]
    ds=0.1,       # Longitudinal loss refinement accuracy [m]
)

loss_loc_refinement.refine_loss_location(result.particles)

# `result.lost_particles` is a snapshot taken right after tracking, so it does
# not see the refined loss locations: re-filter the refined particle object.
# Xcoll flags the different loss mechanisms with distinct non-positive states,
# and the spare capacity kept for secondaries shows up as particle_id < 0.
particles = result.particles
lost_particles = particles.filter((particles.state > -999999999) & (particles.state <= 0))

######################################################
# Plot: Toy ring beam-gas loss map
######################################################
binwidth = 0.1

plt.title(
    f'Toy ring beam-gas loss map '
    f'(beam-gas lifetime: {result.lifetime_tracking/60:.2f} min)'
)
plt.hist(
    lost_particles.s,
    bins=np.arange(0, circumference + binwidth, binwidth),
    weights=lost_particles.weight*1e-3,
)
plt.xlabel('s [m]')
plt.ylabel('Loss rate [kHz]')
plt.grid()
plt.show()