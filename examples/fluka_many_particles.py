import numpy as np
import xpart as xp
import xtrack as xt
import xtrack.particles.pdg as pdg
import xcoll as xc
import time

import matplotlib.pyplot as plt


if xc.FlukaEngine.is_running():
    xc.FlukaEngine.stop()


def run_many_particles(particle_ref, num_part, capacity=None, plot=False):
    if not capacity:
        capacity = num_part*100

    # Create a FLUKA collimator
    coll = xc.FlukaCollimator(length=0.4, assembly='fcc_tcp')
    coll.jaw = 0.001

    # Connect to FLUKA
    xc.FlukaEngine.particle_ref = particle_ref
    xc.FlukaEngine.capacity = capacity
    xc.FlukaEngine.start(elements=coll, clean=True, verbose=False, return_all=True,
                        return_neutral=True, electron_lower_momentum_cut=1.e6)

    # Create an initial distribution of particles, random in 4D, on the left jaw (with the
    # longitudinal coordinates set to zero)
    x_init   = np.random.normal(loc=0.0012, scale=0.2e-3, size=num_part)
    px_init  = np.random.normal(loc=-1.e-5, scale=5.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part_init = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                                   particle_ref=xc.FlukaEngine.particle_ref,
                                   _capacity=xc.FlukaEngine.capacity)
    part = part_init.copy()

    # Do the tracking in FLUKA
    print(f"Tracking {num_part} {pdg.get_name_from_pdg_id(particle_ref.pdg_id[0])}s  (FLUKA)...     ", end='', flush=True)
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
        mass = part.mass[part.pdg_id==pdg_id][0]
        if np.isnan(mass):
            # Need to be careful to avoid NaNs from m0 / m with m = 0
            E = (part.beta0[0]*part.ptau[mask & (part.pdg_id==pdg_id)] + 1)*part.energy0[0]
        else:
            E = part.energy[mask & (part.pdg_id==pdg_id)]
        if len(E[~np.isnan(E)]) > 0:
            en = f"{E[~np.isnan(E)].mean():.1e} ± {E[~np.isnan(E)].std():.1e} eV"
        else:
            en = "NaN"
        print(f"  {num:6} {name:12}{en:21}  (PDG ID: {pdg_id}, mass: {mass} eV)")

    # # Some checks on the FLUKA reference masses for ions. This is useful for the devs.
    # if any([pdg.is_ion(pdg_id) and pdg_id not in xc.fluka_masses for pdg_id in pdg_ids[0]]):
    #     print("New FLUKA reference masses for in database:")
    # from xaux import ProtectFile
    # for pdg_id in pdg_ids[0]:
    #     if pdg.is_ion(pdg_id):
    #         if pdg_id in xc.fluka_masses:
    #             mass = part.mass[part.pdg_id==pdg_id][0]
    #             if not np.isclose(xc.fluka_masses[pdg_id][1], mass):
    #                 print(f"  ERROR: {pdg.get_name_from_pdg_id(pdg_id, long_name=False)}: FLUKA mass in Xcoll {xc.fluka_masses[pdg_id][1]} eV, but got {mass} eV")
    #         else:
    #             mass = part.mass[part.pdg_id==pdg_id][0]
    #             name = f"{pdg.get_name_from_pdg_id(int(pdg_id), long_name=False, subscripts=False)}_MASS_EV_FLUKA"
    #             new_line = f"#define  {name:24}{float(mass):17.4f}    // {int(pdg_id)}   //"
    #             print(new_line)
    #             # Update the reference masses in the database
    #             _ , A, Z, _ = pdg.get_properties_from_pdg_id(pdg_id)
    #             db_file = xc._pkg_root / 'scattering_routines' / 'fluka' / 'reference_masses.py'
    #             with ProtectFile(db_file, 'r+') as fp:
    #                 contents = fp.readlines()
    #                 element_exists_in_table = False
    #                 for i, line in enumerate(contents):
    #                     parts = line.split()
    #                     if line.startswith('#endif /* XCOLL_FLUKA_MASSES_H */'):
    #                         if element_exists_in_table:
    #                             i -= 1
    #                             break
    #                         break
    #                     elif len(parts) > 3 and parts[0] == '#define':
    #                         _, this_A, this_Z, _ = pdg.get_properties_from_pdg_id(parts[4])
    #                         if this_Z == Z:
    #                             element_exists_in_table = True
    #                             if this_A > A:
    #                                 break
    #                         elif this_Z > Z:
    #                             if element_exists_in_table:
    #                                 i -= 1
    #                                 break
    #                             break
    #                 if not element_exists_in_table:
    #                     contents.insert(i, '\n')
    #                 contents.insert(i, new_line + '\n')
    #                 contents = "".join(contents)
    #                 fp.truncate(0)
    #                 fp.seek(0)
    #                 fp.write(contents)
    print()

    # Stop the FLUKA server
    xc.FlukaEngine.stop(clean=True)

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


run_many_particles(xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=6.8e12), 100)
run_many_particles(xt.Particles.reference_from_pdg_id(pdg_id='Pb208', p0c=6.8e12*82), 100, capacity=100_000)
run_many_particles(xt.Particles.reference_from_pdg_id(pdg_id='positron', p0c=200e9), 500, plot=True)


# # Verify reference masses file
# pdg_id = 0
# path = xc._pkg_root / 'scattering_routines' / 'fluka' / 'reference_masses.py'
# with path.open('r') as fid:
#     for line in fid.readlines():
#         if line.startswith('#define') and len(line.split()) > 4:
#             new_pdg_id = int(line.split()[4])
#             assert new_pdg_id > pdg_id
#             # pdg_id = new_pdg_id
#             # mass = pdg.get_mass_from_pdg_id(pdg_id, verbose=False)
#             # print(f"{pdg_id} {(float(line.split()[2]) - mass)/mass*100:.3}%")
