# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xtrack as xt
import xcoll as xc

from many_particles_api import run_many_particles, print_particle_summary, show_new_masses, plot_energy_distribution


def run_many_particles_fluka(particle_ref, num_part, capacity, relative_capacity=None):
    if xc.fluka.engine.is_running():
        xc.fluka.engine.stop()

    # Create a FLUKA collimator
    coll = xc.FlukaCollimator(length=0.1, material='mogr')
    coll.jaw = 0.001

    # Connect to FLUKA
    xc.fluka.engine.particle_ref = particle_ref
    xc.fluka.engine.return_all = True
    xc.fluka.engine.capacity = capacity
    xc.fluka.engine.relative_capacity = relative_capacity
    xc.fluka.engine.start(elements=coll, relative_energy_cut=1e-3, return_all=True, clean=True, verbose=True)

    part = run_many_particles(coll, xc.fluka.engine.particle_ref, num_part, xc.fluka.engine.capacity)

    # Stop the FLUKA server
    xc.fluka.engine.stop(clean=True)

    return part


part = run_many_particles_fluka(xt.Particles('proton', p0c=6.8e12), 1000, capacity=20_000, relative_capacity=100)
print_particle_summary(part)
# show_new_masses(part)

part = run_many_particles_fluka(xt.Particles('Pb208', p0c=6.8e12*82), 50, capacity=100_000, relative_capacity=2000)
print_particle_summary(part)
# show_new_masses(part)

part = run_many_particles_fluka(xt.Particles('positron', p0c=200e9), 1000, capacity=50_000, relative_capacity=100)
print_particle_summary(part)
# show_new_masses(part)
plot_energy_distribution(part, nbins=250)
