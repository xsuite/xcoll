# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xtrack as xt
import xcoll as xc

from many_particles_api import run_many_particles, print_particle_summary, show_new_masses, plot_energy_distribution


def run_many_particles_geant4(particle_ref, num_part, capacity):
    if xc.geant4.engine.is_running():
        xc.geant4.engine.stop()

    # Create a Geant4 collimator
    coll = xc.Geant4Collimator(length=0.1, material='mogr')
    coll.jaw = 0.001

    # Connect to Geant4
    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.return_all = True
    xc.geant4.engine.start(elements=coll, relative_energy_cut=1e-3, return_all=True, clean=True, verbose=True)

    part = run_many_particles(coll, xc.geant4.engine.particle_ref, num_part, capacity)

    # Stop the Geant4 server
    xc.geant4.engine.stop(clean=True)

    return part


part = run_many_particles_geant4(xt.Particles('proton',   p0c=6.8e12),    25_000, capacity=500_000)
print_particle_summary(part)
# show_new_masses(part)

part = run_many_particles_geant4(xt.Particles('Pb208',    p0c=6.8e12*82),   2000, capacity=500_000)
print_particle_summary(part)
# show_new_masses(part)

part = run_many_particles_geant4(xt.Particles('positron', p0c=200e9),     25_000, capacity=500_000)
print_particle_summary(part)
# show_new_masses(part)
plot_energy_distribution(part, nbins=250)
