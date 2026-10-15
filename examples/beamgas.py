# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

import xobjects as xo
import xtrack as xt

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
# Install beam-gas scattering markers
######################################################
# Each marker represents the section of the line from the previous marker to
# itself: they are installed at equally spaced positions, the last one at the
# end of the line so that the sections cover the whole ring.
circumference = line.get_length()
n_beamgas = 16

line.xcoll.beamgas.install_markers(
    at=list(np.linspace(0, circumference, n_beamgas + 1)[1:]))

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
tt_beamgas = tab.rows[tab.element_type == 'BeamGasMarker']
beamgas_markers = tt_beamgas.name

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
# Configure the beam-gas markers
######################################################
line.discard_tracker()
line.build_tracker(_context=xo.ContextCpu(omp_num_threads='auto'))

# Integrates the gas density over the section of each marker, and stores on
# the marker the interaction rate together with the local gas composition,
# optics, beam parameters and Monte Carlo settings. The configured line can
# be saved here (e.g. with line.to_json): once loaded again, it generates the
# same particles, e.g. in separate jobs.
line.xcoll.beamgas.configure(
    gas_density=gas_density,
    process='coulomb',
    # Losses from angles outside this window are missing from the result, so
    # it must extend well below the smallest angle that causes a loss (here
    # ~3.5 mrad) and well above the angle beyond which every particle is lost
    # (here ~20 mrad). The angles are sampled log-uniformly with importance
    # weights, so a wide window only costs statistics as
    # log(theta_max/theta_min). BeamGasResult.cutoff_scan and
    # rate_above_theta_max show whether the window is wide enough.
    coulomb_theta=(1e-6, 0.3),
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

print(line.xcoll.beamgas.analyse().local_rates)

######################################################
# Generate and track the beam-gas-scattered particles
######################################################
# For each marker:
#   1. Generate the particles scattered in its section, each with one
#      interaction with the residual gas, and the record of the interactions
#   2. Track them around the ring for `nturns` turns, starting and ending at
#      the marker
nturns = 100
particles_by_marker = {}
logs = {}
for name in beamgas_markers:
    particles, log = line.xcoll.beamgas.generate(at=name, return_log=True)
    line.track(particles, ele_start=name, ele_stop=name, num_turns=nturns,
               with_progress=1)
    particles_by_marker[name] = particles
    logs[name] = log

######################################################
# Beam-gas loss rate and lifetime
######################################################
# Marker by marker, with the records of the interactions: the result also
# contains the Monte Carlo error and the scan of the lower generation cut.
# The merged particles alone give the loss rate and the lifetime.
result = line.xcoll.beamgas.analyse(particles_by_marker, logs=logs)

# All of the above in one call:
# result = line.xcoll.beamgas.run(n_turns=nturns)

print(f'Beam-gas interaction rate: {result.rate_scattering*1e-3:.3f} kHz')
print(f'Beam-gas loss rate:        {result.rate_tracking*1e-3:.3f} kHz')
print(f'Beam-gas lifetime:         {result.lifetime_tracking/60:.2f} '
      f'+- {result.lifetime_tracking_error/60:.2f} min')

# Lifetime that a lower cut `cut` would have given: it must be flat in the
# first rows, otherwise coulomb_theta[0] is too high
print(result.cutoff_scan)
print(f'Interaction rate above theta_max: '
      f'{result.rate_above_theta_max/result.rate_tracking:.1e} of the loss rate')

# The weights are rates [1/s]: merged, the samples of all the markers
# represent the beam-gas scattering of the whole ring
particles = xt.Particles.merge(list(particles_by_marker.values()))

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

loss_loc_refinement.refine_loss_location(particles)

# Xcoll flags the different loss mechanisms with distinct non-positive states,
# and the spare capacity kept for secondaries shows up as particle_id < 0.
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