import numpy as np
import xpart as xp
import xtrack as xt
import xtrack.particles.pdg as pdg
import xcoll as xc
from  xcoll import constants as xcs
import time

import matplotlib.pyplot as plt

# # Treat warnings as errors to debug
# np.set_printoptions(threshold=np.inf)
# import warnings
# warnings.filterwarnings("error")


if xc.geant4.engine.is_running():
    xc.geant4.engine.stop()


def run_many_particles(particle_ref, num_part, capacity=None, plot=False):

    # Create a Geant4 collimator
    coll = xc.Geant4Collimator(length=0.4, material='mogr')
    coll.jaw = 0.001

    # Connect to Geant4
    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=coll, clean=True, verbose=False)

    # Create an initial distribution of particles, random in 4D, on the left jaw (with the
    # longitudinal coordinates set to zero)
    x_init   = np.random.normal(loc=0.002, scale=0.2e-3, size=num_part)
    px_init  = np.random.normal(loc=-1.e-5, scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                   particle_ref=xc.geant4.engine.particle_ref,
                                   _capacity=capacity)
    part = part_init.copy()

    # Do the tracking in Geant4
    print(f"Tracking {num_part} {pdg.get_name_from_pdg_id(particle_ref.pdg_id[0])}s...     ", end='', flush=True)
    start = time.time()
    coll.track(part)
    print(f"Done in {round(time.time()-start, 3)}s.", flush=True)

    # Get only the initial particles that survived and all new particles (even if dead, as neutral particles will be flagged dead)
    mask = (part.state > 0) | (part.particle_id >= num_part)
    pdg_ids = np.unique(part.pdg_id[mask], return_counts=True)
    print("Returned particles:")
    for pdg_id, num in zip(*pdg_ids):
        try:
            name = pdg.get_name_from_pdg_id(pdg_id, long_name=False)
        except ValueError:
            name = 'unknown'
        if part.state[part.pdg_id==pdg_id][0] == xcs.MASSLESS_OR_NEUTRAL:
            mass = 0
        else:
            mass = part.mass[part.pdg_id==pdg_id][0]
        E = part.energy[mask & (part.pdg_id==pdg_id)]
        en = f"{E[~np.isnan(E)].mean():.1e} ± {E[~np.isnan(E)].std():.1e} eV"
        print(f"  {num:6} {name:12}{en:21}  (PDG ID: {pdg_id}, mass: {mass} eV)")
    print()

    # Stop the Geant4 server
    xc.geant4.engine.stop(clean=True)

    if plot:
        data = []
        minimum = 1e10
        maximum = 0
        for pdg_id in pdg_ids[0]:
            if np.isnan(part.mass[part.pdg_id==pdg_id][0]):
                # Need to be careful to avoid NaNs from m0 / m with m = 0
                data.append((part.beta0[0]*part.ptau[part.pdg_id == pdg_id] + 1)*part.energy0[0])
            else:
                data.append(part.energy[part.pdg_id == pdg_id])
            minimum = min(minimum, np.min(data[-1]))
            maximum = max(maximum, np.max(data[-1]))
        bins = np.logspace(np.log10(minimum), np.log10(maximum), 50)
        plt.figure(figsize=(8, 6))
        for pdg_id, this_data in zip(pdg_ids[0], data):
            try:
                name = pdg.get_name_from_pdg_id(pdg_id, long_name=True)
            except ValueError:
                name = 'unknown'
            plt.hist(this_data, bins=bins, density=True, alpha=0.4, label=name)

        # Set horizontal axis to logarithmic scale
        plt.xscale('log')
        plt.yscale('log')

        # Labels and legend
        plt.xlabel('Energy [eV]')
        plt.ylabel('Normalised frequency')
        plt.legend()
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.show()


run_many_particles(xt.Particles('proton', p0c=6.8e12),   10_000, capacity=20_000)
run_many_particles(xt.Particles('Pb208', p0c=6.8e12*82), 5_000, capacity=10_000)
run_many_particles(xt.Particles('positron', p0c=200e9),  10_000, capacity=20_000, plot=True)
