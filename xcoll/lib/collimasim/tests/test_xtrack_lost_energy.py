import numpy as np
import xtrack as xt
import xpart as xp

import collimasim as cs
import pandas as pd

def test_energy_diffrential():
    N_part = 10
    N_max_expected_products = 10

    particles = xp.Particles(
        _capaccity = N_part,
		p0c = 36900e9,
        x = np.zeros(N_part),
		y = np.zeros(N_part),
		zeta = np.zeros(N_part),
		px = np.zeros(N_part),
		py = np.zeros(N_part),
		delta = np.zeros(N_part),
        pdg_id = xp.pdg.get_pdg_id_from_name('Pb208'),
        mass0=193687676197.13638,
        )

    particles.x[1] = 0.04226982005192951 + 10e-6
    particles.px[1] = 1e-7
    particles.x[2] = -0.04226982005192951 - 10e-6
    particles.px[2] = -1e-7
    particles.x[3] = 0.04226982005192951 + 10e-6
    particles.y[4] = 0.05

    part_orig = particles.copy()

    # Random parameters to get some scattering in the collimators
    g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_forions.dat",
                                        bdsim_config_file="resources/settings_ions.gmad",
                                        tfs_file="resources/collonly_twiss_file_example.tfs",
                                        reference_pdg_id=xp.pdg.get_pdg_id_from_name('Pb208'),
                                        reference_kinetic_energy=36706.8206525e9,
                                        emittance_norm=(0.000521429683675842, 0.000521429683675842),
                                        relative_energy_cut=0.0001,
                                        seed=1993,
                                        batchMode=True)

    coll1 = cs.Geant4Collimator(name="tcp.a.b1", g4manager=g4man)


    prods = coll1.interact(particles)

    part_df = particles.to_pandas(compact=True)
    part_orig_df = part_orig.to_pandas(compact=True)

    part_mass_ratio = part_df.charge_ratio / part_df.chi
    part_mom = (part_df.delta + 1) * part_df.p0c * part_mass_ratio
    part_mass = part_mass_ratio * part_df.mass0
    part_tot_energy = np.sqrt(part_mom**2 + part_mass**2)
    part_df['energy'] = part_tot_energy

    part_orig_mass_ratio = part_orig_df.charge_ratio / part_orig_df.chi
    part_orig_mom = (part_orig_df.delta + 1) * part_orig_df.p0c * part_orig_mass_ratio
    part_orig_mass = part_orig_mass_ratio * part_orig_df.mass0
    part_orig_tot_energy = np.sqrt(part_orig_mom**2 + part_orig_mass**2)
    part_orig_df['energy'] = part_tot_energy

    prod_df = pd.DataFrame(prods)

    prod_mass_ratio = prod_df.mass_ratio
    prod_mom = (prod_df.delta + 1) * part_df.p0c[0] * prod_mass_ratio
    prod_mass = prod_mass_ratio * part_df.mass0[0]
    prod_tot_energy = np.sqrt(prod_mom**2 + prod_mass**2)
    prod_df['energy'] = prod_tot_energy

    for part_id in [1, 2, 3]:
        pp = part_df[part_df['particle_id'] == part_id]
        pp_orig = part_df[part_orig_df['particle_id'] == part_id]
        prods_pp = prod_df[prod_df['parent_particle_id'] == part_id]

        combined_energy = float(pp['energy']) + float(sum(prods_pp['energy']))
        energy_diff = combined_energy - float(pp_orig['energy'])



def main():
    test_energy_diffrential()


if __name__ == '__main__':
    main()
