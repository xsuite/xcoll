# This is a script to test the FLUKA integration into xtrack.
# Before running the script, the source code needs to be compiled. In the root
# package folder, run ./compile_fluka.sh to do so.

import numpy as np
import xpart as xp
import xtrack as xt
import xcoll as xc
import time

import matplotlib.pyplot as plt


num_part = 100_000
capacity = 2*num_part
particle_ref = xt.Particles('proton', p0c=450e9)
material = xc.materials.MolybdenumGraphite
length = 0.2

impact_parameter = 2e-3
impact_spread = 0.2e-6
angular_divergence = 5.e-6
nbins = 200
E_min = 50e9
E_zoom = 350e9
E_high = particle_ref.energy0[0]


# Stop any running engines
if xc.fluka.engine.is_running():
    xc.fluka.engine.stop(clean=True)
if xc.geant4.engine.is_running():
    xc.geant4.engine.stop(clean=True)


# Create collimators in Everest, FLUKA, and Geant4
jaw = 0.001
coll1 = xc.EverestCollimator(length=length, material=material, jaw=jaw)
coll2 = xc.FlukaCollimator(length=length, material=material, jaw=jaw)
coll3 = xc.Geant4Collimator(length=length, material=material, jaw=jaw)


# Connect to FLUKA
xc.fluka.engine.particle_ref = particle_ref
xc.fluka.engine.capacity = capacity
xc.fluka.engine.seed = 5656565
xc.fluka.engine.return_none = True
xc.fluka.engine.return_protons = True
xc.fluka.engine.start(elements=coll2, clean=True, verbose=False)
xc.fluka.engine.init_tracking(max_particle_id=num_part)  # Prepares FLUKA for tracking, needed because we will track in chunks

# Connect to Geant4
xc.geant4.engine.particle_ref = particle_ref
xc.geant4.engine.seed = 5656565
xc.geant4.engine.return_none = True
xc.geant4.engine.return_protons = True
xc.geant4.engine.start(elements=coll3, clean=True, verbose=False)


# Create an initial distribution of particles, random in 4D, on the left jaw (with the
# longitudinal coordinates set to zero)
np.random.seed(seed=23823842)
x_init   = np.random.normal(loc=jaw+impact_parameter, scale=impact_spread, size=num_part)
px_init  = np.random.normal(loc=0., scale=angular_divergence, size=num_part)
y_init   = np.random.normal(loc=0., scale=impact_spread, size=num_part)
py_init  = np.random.normal(loc=0., scale=angular_divergence, size=num_part)
part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                               particle_ref=xc.fluka.engine.particle_ref,
                               _capacity=capacity)


# Do the tracking in Everest
start = time.time()
energy_everest = np.array([], dtype=np.float64)
# We split it in small chunks to be more efficient
num_step = min(1_000_000, num_part)
steps = num_part // num_step
for i in range(steps):
    print(f"Tracking {num_part} particles (Everest)...  Step {i}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[i*num_step:(i+1)*num_step] = True
    part = part_init.filter(mask)
    coll1.track(part)
    this_energy_everest = part.energy[part.state > 0]
    energy_everest = np.concatenate((energy_everest, this_energy_everest))
if num_part > steps*num_step:  # Some leftover particles
    print(f"Tracking {num_part} particles (Everest)...  Step {steps}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[steps*num_step:num_part] = True
    part = part_init.filter(mask)
    coll1.track(part)
    this_energy_everest = part.energy[part.state > 0]
    energy_everest = np.concatenate((energy_everest, this_energy_everest))
print(f"Tracking {num_part} particles (Everest)...  Done in {round(1000*(time.time()-start), 3)}ms")


# Do the tracking in FLUKA
start = time.time()
energy_fluka = np.array([], dtype=np.float64)
# We split it in small chunks to be more efficient
num_step = min(10_000, num_part)
steps = num_part // num_step
for i in range(steps):
    print(f"Tracking {num_part} particles (FLUKA)...    Step {i}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[i*num_step:(i+1)*num_step] = True
    mask[num_part + i*num_step:num_part + (i+1)*num_step] = True  # For secondaries
    part = part_init.filter(mask)
    coll2.track(part)
    mask = (part.state > 0)
    parents, inv = np.unique(part.parent_particle_id[mask], return_inverse=True)
    this_energy_fluka = np.zeros_like(parents, dtype=np.float64)
    np.add.at(this_energy_fluka, inv, part.energy[mask])
    energy_fluka = np.concatenate((energy_fluka, this_energy_fluka))
if num_part > steps*num_step:  # Some leftover particles
    print(f"Tracking {num_part} particles (FLUKA)...    Step {steps}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[steps*num_step:num_part] = True
    mask[num_part + steps*num_step:num_part + num_part] = True  # For secondaries
    part = part_init.filter(mask)
    coll2.track(part)
    mask = (part.state > 0)
    parents, inv = np.unique(part.parent_particle_id[mask], return_inverse=True)
    this_energy = np.zeros_like(parents, dtype=np.float64)
    np.add.at(this_energy, inv, part.energy[mask])
    energy_fluka = np.concatenate((energy_fluka, this_energy))
print(f"Tracking {num_part} particles (FLUKA)...    Done in {round(time.time()-start, 3)}s.")


# Do the tracking in Geant4
start = time.time()
energy_geant4 = np.array([], dtype=np.float64)
# We split it in small chunks to be more efficient
num_step = min(10_000, num_part)
steps = num_part // num_step
for i in range(steps):
    print(f"Tracking {num_part} particles (Geant4)...   Step {i}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[i*num_step:(i+1)*num_step] = True
    mask[num_part + i*num_step:num_part + (i+1)*num_step] = True  # For secondaries
    part = part_init.filter(mask)
    coll3.track(part)
    mask = (part.state > 0)
    parents, inv = np.unique(part.parent_particle_id[mask], return_inverse=True)
    this_energy = np.zeros_like(parents, dtype=np.float64)
    np.add.at(this_energy, inv, part.energy[mask])
    energy_geant4 = np.concatenate((energy_geant4, this_energy))
if num_part > steps*num_step:  # Some leftover particles
    print(f"Tracking {num_part} particles (Geant4)...   Step {steps}/{steps}", flush=True, end='\r')
    mask = np.zeros(capacity, dtype=np.bool)
    mask[steps*num_step:num_part] = True
    mask[num_part + steps*num_step:num_part + num_part] = True  # For secondaries
    part = part_init.filter(mask)
    coll3.track(part)
    mask = (part.state > 0)
    parents, inv = np.unique(part.parent_particle_id[mask], return_inverse=True)
    this_energy = np.zeros_like(parents, dtype=np.float64)
    np.add.at(this_energy, inv, part.energy[mask])
    energy_geant4 = np.concatenate((energy_geant4, this_energy))
print(f"Tracking {num_part} particles (Geant4)...   Done in {round(time.time()-start, 3)}s.")


# Stop the engines
xc.fluka.engine.stop(clean=True)
xc.geant4.engine.stop(clean=True)


print(f"Deposited energy (Everest): {E_high - energy_everest.sum()/num_part:.3e} eV per particle")
print(f"Deposited energy (FLUKA):   {E_high - energy_fluka.sum()/num_part:.3e} eV per particle")
print(f"Deposited energy (Geant4):  {E_high - energy_geant4.sum()/num_part:.3e} eV per particle")


# Create bins based on the min and max energy
bins = np.logspace(np.log10(E_min), np.log10(E_high), nbins + 1)
bin_centres = np.sqrt(bins[:-1] * bins[1:])
log_edges = np.log10(bins)
dlog = np.diff(log_edges)
bins_zoom = np.logspace(np.log10(E_zoom), np.log10(E_high), nbins + 1)
bin_centres_zoom = np.sqrt(bins_zoom[:-1] * bins_zoom[1:])
log_edges_zoom = np.log10(bins_zoom)
dlog_zoom = np.diff(log_edges_zoom)

fig, ax = plt.subplots(2, 1, figsize=(14, 7))
for ee, label in zip([energy_everest, energy_fluka, energy_geant4],
                     ['Everest', 'FLUKA', 'Geant4']):
    # Full plot
    counts, edges = np.histogram(ee, bins=bins)
    dNdlogE = counts / (len(ee) * dlog)
    ax[0].step(bin_centres, dNdlogE, where='mid', label=label)
    max_dNdlogE = dNdlogE.max()
    # Zoom plot
    counts, edges = np.histogram(ee, bins=bins_zoom)
    dNdlogE = counts / (len(ee) * dlog_zoom)
    ax[1].step(bin_centres_zoom, dNdlogE, where='mid', label=label)

# Set horizontal axis to logarithmic scale
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_yscale('log')

# Labels and legend
ax[0].set_xlabel('Energy [eV]')
ax[0].set_ylabel('Normalised frequency ' + r'$\frac{dN}{d\log E}$')
ax[0].legend(loc='upper left')
ax[0].grid(True, which='both', linestyle='--', linewidth=0.5)

ax[1].set_xlabel('Energy [eV]')
ax[1].set_ylabel('Normalised frequency ' + r'$\frac{dN}{d\log E}$')
ax[1].legend(loc='upper left')
ax[1].grid(True, which='both', linestyle='--', linewidth=0.5)

text = f"""\
Deposited {E_high - energy_everest.sum()/num_part:.3e} eV per particle (Everest)
Deposited {E_high - energy_fluka.sum()/num_part:.3e} eV per particle (FLUKA)
Deposited {E_high - energy_geant4.sum()/num_part:.3e} eV per particle (Geant4)"""
print(max_dNdlogE)
ax[0].text(E_min*(E_high/E_min)**0.07, max_dNdlogE, text, va='top',
           fontsize=10, bbox=dict(edgecolor='k', facecolor='w', alpha=0.5, boxstyle='round'))

plt.suptitle(f"Scattering {num_part} protons at {E_high:.3e} eV on {length}m {material.name}")
plt.tight_layout()
plt.savefig('plots/compare_energy_3codes.png', dpi=300)
plt.show()
