import numpy as np
import xpart as xp
import xtrack as xt
import xtrack.particles.pdg as pdg
import xcoll as xc
from  xcoll import constants as xcs
import time

import matplotlib.pyplot as plt

from xcoll.xoconstants import constant

# # Treat warnings as errors to debug
# np.set_printoptions(threshold=np.inf)
# import warnings
# warnings.filterwarnings("error")


if xc.geant4.engine.is_running():
    xc.geant4.engine.stop()


def run_many_particles(particle_ref, num_part, capacity=None, plot=False):

    # Create a Geant4 collimator
    coll = xc.Geant4Collimator(length=0.1, material='mogr')
    coll.jaw = 0.001

    # Connect to Geant4
    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=coll, relative_energy_cut=1e-3, return_all=True, clean=True, verbose=True)

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

    # Stop the Geant4 server
    xc.geant4.engine.stop(clean=True)

    return part


def print_particle_summary(part):
    num_part = part.particle_id[part.particle_id == part.parent_particle_id].max()
    # Get only the initial particles that survived and all new particles (even if dead, as neutral particles will be flagged dead)
    mask = (part.state > 0) | (part.particle_id >= num_part)
    pdg_ids = np.unique(part.pdg_id[mask], return_counts=True)
    idx = np.argsort(pdg_ids[1])[::-1]
    print(f"Returned {np.sum(mask)} particles ({np.sum(part.particle_id >= num_part)} secondaries):")
    for pdg_id, num in zip(pdg_ids[0][idx], pdg_ids[1][idx]):
        try:
            name = pdg.get_name_from_pdg_id(pdg_id, long_name=False)
        except ValueError:
            name = 'unknown'
        if part.state[part.pdg_id==pdg_id][0] == xcs.MASSLESS_OR_NEUTRAL:
            mass = 0
        else:
            mass = part.mass[part.pdg_id==pdg_id][0]
        E = part.energy[mask & (part.pdg_id==pdg_id)]*1e-9
        q1 = np.percentile(E[~np.isnan(E)], 25)
        med = np.percentile(E[~np.isnan(E)], 50)
        q3 = np.percentile(E[~np.isnan(E)], 75)
        en = f"{med:.3e} ∊ [{q1:.1e}, {q3:.1e}] GeV  (50% ± 25%)"
        print(f"  {num:6} {name:12}{en:21}    (PDG ID: {pdg_id:5}, mass: {mass*1e-9:.4f} GeV)")
    print()


def plot_energy_distribution(part, *, nbins=50, max_types=6):
    num_part = part.particle_id[part.particle_id == part.parent_particle_id].max()
    # Get only the initial particles that survived and all new particles (even if dead, as neutral particles will be flagged dead)
    mask = (part.state > 0) | (part.particle_id >= num_part)
    pdg_ids = np.unique(part.pdg_id[mask], return_counts=True)
    idx = np.argsort(pdg_ids[1])[::-1][:max_types]
    data = []
    minimum = 1e10
    maximum = 0
    for pdg_id in pdg_ids[0][idx]:
        if np.isnan(part.mass[part.pdg_id==pdg_id][0]):
            # Need to be careful to avoid NaNs from m0 / m with m = 0
            data.append((part.beta0[0]*part.ptau[part.pdg_id == pdg_id] + 1)*part.energy0[0])
        else:
            data.append(part.energy[part.pdg_id == pdg_id])
        minimum = min(minimum, np.min(data[-1]))
        maximum = max(maximum, np.max(data[-1]))
    bins = np.logspace(np.log10(minimum), np.log10(maximum), nbins + 1)
    bin_centres = np.sqrt(bins[:-1] * bins[1:])
    N_tot = sum(len(arr) for arr in data)
    log_edges = np.log10(bins)
    dlog = np.diff(log_edges)

    # plt.figure(figsize=(8, 6))
    fig, ax = plt.subplots(2, 1, figsize=(14, 7))
    for pdg_id, this_data in zip(pdg_ids[0][idx], data):
        try:
            name = pdg.get_name_from_pdg_id(pdg_id, long_name=True)
        except ValueError:
            name = 'unknown'
        # ax[0].hist(this_data, bins=bins, density=True, alpha=0.4, label=name)
        counts, edges = np.histogram(this_data, bins=bins)
        bin_widths = np.diff(edges)
        dNdlogE = counts / (N_tot * dlog) # Global normalisation: divide by N_tot, not by len(E)
        ax[0].step(bin_centres, dNdlogE, where='mid', label=name)

        dNdE = counts / (N_tot * bin_widths)  # Global normalisation: divide by N_tot, not by len(E)
        ax[1].step(bin_centres, dNdE, where='mid', label=name)

    # Set horizontal axis to logarithmic scale
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')

    # Labels and legend
    ax[0].set_xlabel('Energy [eV]')
    ax[0].set_ylabel('Normalised frequency ' + r'$\frac{dN}{d\log E}$')
    ax[0].legend()
    ax[1].set_xlabel('Energy [eV]')
    ax[1].set_ylabel('Normalised frequency ' + r'$\frac{dN}{dE}$')
    ax[1].legend()
    ax[0].grid(True, which='both', linestyle='--', linewidth=0.5)
    ax[1].grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()


def show_new_masses(part, show_existing=False):
    mask = part.particle_id != part.parent_particle_id
    pdg_ids = np.unique(part.pdg_id[mask])
    mess_new = []
    mess_exist = []
    for pdg_id in pdg_ids:
        if pdg_id == -999999999:
            continue
        if part.state[part.pdg_id==pdg_id][0] == xcs.MASSLESS_OR_NEUTRAL:
            continue
        try:
            if pdg.is_ion(pdg_id):
                name = pdg.get_name_from_pdg_id(pdg_id, long_name=False, subscripts=False)
            else:
                name = pdg.get_name_from_pdg_id(pdg_id, long_name=True).upper()
        except ValueError:
            name = 'unknown'
        masses = part.mass[mask & (part.pdg_id==pdg_id)]
        q1  = np.percentile(masses[~np.isnan(masses)], 25)
        med = np.percentile(masses[~np.isnan(masses)], 50)
        q3  = np.percentile(masses[~np.isnan(masses)], 75)
        if abs(pdg_id) not in xc.geant4.particle_masses:
            mess  = f"New mass for PDG ID {pdg_id} ({name}):  {med} eV "
            mess += f"[{masses.min()-med:.2e}, {q1-med:.2e}, +{q3-med:.2e}, +{masses.max()-med:.2e}] "
            mess += f"({len(masses)} samples)."
            mess_new.append(mess)
        else:
            mass  = xc.geant4.particle_masses[pdg_id]
            mess  = f"Existing mass for PDG ID {pdg_id} ({name}):  {mass} eV (existing) vs {med} eV (returned). "
            mess += f"Difference {med - mass:.2e} eV "
            mess += f"({len(masses)} samples)."
            mess_exist.append(mess)
    for mess in mess_new:
        print(mess)
    if show_existing:
        for mess in mess_exist:
            print(mess)


part = run_many_particles(xt.Particles('proton',   p0c=6.8e12),    25_000, capacity=500_000)
print_particle_summary(part)
# show_new_masses(part)

part = run_many_particles(xt.Particles('Pb208',    p0c=6.8e12*82),   2000, capacity=500_000)
print_particle_summary(part)
# show_new_masses(part)

part = run_many_particles(xt.Particles('positron', p0c=200e9),     25_000, capacity=500_000)
print_particle_summary(part)
# show_new_masses(part)
plot_energy_distribution(part, nbins=250)
