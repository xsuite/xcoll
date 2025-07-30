import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp

import collimasim as cs
from matplotlib import pyplot as plt

import multiprocessing

def make_particles():
    np.random.seed(seed=1994)
    n_part = 100
    particles = xp.Particles(
        _capacity=n_part,
		p0c = 1.e3, # 1 KeV dummy particles
        mass0 = xp.ELECTRON_MASS_EV,
        x=np.linspace(0.1, 0.3, n_part), # Offsets are x=0.2, y=0.3
        #x=np.random.uniform(0.1, 0.3, n_part), # Offsets are x=0.2, y=0.3
        px=np.zeros(n_part),
        y=np.linspace(0.2, 0.4, n_part),
        #y=np.random.uniform(0.2, 0.4, n_part),
        #y=np.zeros(n_part) + 0.3,
        py=np.zeros(n_part),
        zeta=np.zeros(n_part),
        delta=np.zeros(n_part),
    )

    return particles


def test_xtrack_tilt():
    np.random.seed(seed=1994)
    particles = make_particles()

    g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_testing.dat",
                                        bdsim_config_file="resources/settings_black_absorber.gmad",
                                        tfs_file="resources/twiss_file_testing.tfs",
                                        reference_pdg_id=-11,
                                        reference_kinetic_energy=0.510998946e6 * 10000,
                                        emittance_norm=(1.e-6, 1.e-6),
                                        relative_energy_cut=0.001,
                                        seed=1993,
                                        # batchMode=False
                                        batchMode=True
                                        )
    
    coll1 = g4man.make_xtg4_collimator("coll_skew_tilted") # Use new convenience method

    # Generate a simple sequence
    line = xt.Line(
        elements=[coll1,
                  ])

    context = xo.ContextCpu()
    line.build_tracker(_context=context)

    part_copy = particles.copy()
    line.track(part_copy, num_turns=1)
    
    print(f"Lost particles: {sum(part_copy.state==-333)} / {len(part_copy.state)}")
    assert sum(part_copy.state==-333) == 89

    return part_copy


def main():
    part_g4 = test_xtrack_tilt()

    # plt.scatter(part_g4.x, part_g4.y, c=part_g4.state)
    # ax = plt.gca()
    # ax.axvline(0.15, c='r')
    # ax.axvline(0.25, c='r')
    # plt.show()

    print('Done!')


if __name__ == '__main__':
    main()
