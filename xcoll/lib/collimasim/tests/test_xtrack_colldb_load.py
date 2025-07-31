import numpy as np
import collimasim as cs

def test_yaml_colldb_load():
    try: # yaml file loading is not implemented yet, so expect an error
        g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_yaml_example.yaml",
                                            bdsim_config_file="resources/settings_black_absorber.gmad",
                                            tfs_file="resources/twiss_file_testing.tfs",
                                            reference_pdg_id=-11,
                                            reference_kinetic_energy=182.4994890018054e9, # eV
                                            emittance_norm=(0.000521429683675842, 1.0357164949725628e-06),
                                            relative_energy_cut=0.001,
                                            seed=1993,
                                            batchMode=True)
    except ValueError as e:
        print(f"YAML file loading sucessfuly reports an error: {e}")


def test_old_colldb_load():
    g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_old_example.dat",
                                        bdsim_config_file="resources/settings_black_absorber.gmad",
                                        tfs_file="resources/collonly_twiss_file_example.tfs",
                                        reference_pdg_id=-11,
                                        reference_kinetic_energy=182.4994890018054e9, # eV
                                        emittance_norm=(0.000521429683675842, 1.0357164949725628e-06),
                                        relative_energy_cut=0.001,
                                        seed=1993,
                                        batchMode=True)
    assert np.isclose(g4man.collimators['tcp.a.b1']['halfgap'], 0.00976272689045678, atol=1e-10, rtol=1e-6)


def test_new_colldb_load():
    g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_new_example.dat",
                                        bdsim_config_file="resources/settings_black_absorber.gmad",
                                        tfs_file="resources/collonly_twiss_file_example.tfs",
                                        reference_pdg_id=-11,
                                        reference_kinetic_energy=182.4994890018054e9, # eV
                                        emittance_norm=(0.000521429683675842, 1.0357164949725628e-06),
                                        relative_energy_cut=0.001,
                                        seed=1993,
                                        batchMode=True)
    assert np.isclose(g4man.collimators['tcp.a.b1']['halfgap'], 0.00976272689045678, atol=1e-10, rtol=1e-6)


def test_collgaps_load():
    g4man = cs.Geant4CollimationManager(collimator_file="resources/collgaps.dat",
                                        bdsim_config_file="resources/settings_black_absorber.gmad",
                                        tfs_file="resources/collonly_twiss_file_example.tfs",
                                        reference_pdg_id=-11,
                                        reference_kinetic_energy=182.4994890018054e9, # eV
                                        emittance_norm=(0.000521429683675842, 1.0357164949725628e-06),
                                        relative_energy_cut=0.001,
                                        seed=1993,
                                        batchMode=True)
    assert np.isclose(g4man.collimators['tcp.a.b1']['halfgap'], 0.00976272689045678, atol=1e-10, rtol=1e-6)


def test_badfile_load():
    try:
        # This is an example of loading an invalid file
        g4man = cs.Geant4CollimationManager(collimator_file="resources/settings.gmad",
                                            bdsim_config_file="resources/settings_black_absorber.gmad",
                                            tfs_file="resources/collonly_twiss_file_example.tfs",
                                            reference_pdg_id=-11,
                                            reference_kinetic_energy=182.4994890018054e9, # eV
                                            emittance_norm=(0.000521429683675842, 1.0357164949725628e-06),
                                            relative_energy_cut=0.001,
                                            seed=1993,
                                            batchMode=True)
    except ValueError as e:
         print(f"Invalid file loading sucessfuly reports an error: {e}")


def main():
    # test_yaml_colldb_load()
    test_new_colldb_load()
    # test_old_colldb_load()
    # test_collgaps_load()
    # test_badfile_load()

if __name__ == '__main__':
    main()


