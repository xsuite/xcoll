import numpy as np
import xtrack as xt
import xpart as xp

import collimasim as cs


def test_object_creation():
    N_part = 10
    N_max_expected_products = 10

    tot_num_part = N_part

    particles = xp.Particles(
		p0c = 36900e9,
        x = np.zeros(tot_num_part),
		y = np.zeros(tot_num_part),
		zeta = np.zeros(tot_num_part),
		px = np.zeros(tot_num_part),
		py = np.zeros(tot_num_part),
		delta = np.zeros(tot_num_part))

    particles.num_particles = N_part

    particles.x[1] = 0.3
    particles.y[2] = 0.3

    # Random parameters to get some scattering in the collimators
    g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_forions.dat",
                                        bdsim_config_file="resources/settings_ions.gmad",
                                        tfs_file="resources/collonly_twiss_file_example.tfs",
                                        reference_pdg_id=xp.pdg.get_pdg_id_from_name('Pb208'),
                                        reference_kinetic_energy=36706.8206525e9,
                                        emittance_norm=(0.000521429683675842, 0.000521429683675842),
                                        relative_energy_cut=0.001,
                                        seed=1993,
                                        batchMode=True)

    coll1 = cs.Geant4Collimator(name="tcp.a.b1", g4manager=g4man)
    coll2 = cs.Geant4Collimator(name="tcp.b.b1", g4manager=g4man)


    print("Particle coordinates (x) before first collimator")
    print (particles.delta)

    prods = coll1.interact(particles)

    print("Particle coordinates (x) after first collimator")
    print (particles.x)
    print (particles.y)

    print("="*30)
    print("Products:")
    for crd in prods:
        print(crd, ":", prods[crd])
    print("="*30)

    prods = coll2.interact(particles)

    print("Particle coordinates (x) after second collimator")
    print (particles.x)
    print (particles.y)

    print("="*30)
    print("Products:")
    for crd in prods:
        print(crd, ":", prods[crd])
    print("="*30)

    assert len(prods['x']) == 59


def main():
    test_object_creation()


if __name__ == '__main__':
    main()
