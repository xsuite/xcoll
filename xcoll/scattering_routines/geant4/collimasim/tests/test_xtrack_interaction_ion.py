import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp

import collimasim as cs

import multiprocessing

def test_fodo():
    np.random.seed(seed=1994)

    # Random parameters to get some scattering in the collimators
    g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_forions.dat",
                                        bdsim_config_file="resources/settings_ions.gmad",
                                        tfs_file="resources/collonly_twiss_file_example.tfs",
                                        reference_pdg_id=xp.pdg.get_pdg_id_from_name('Pb208'),
                                        reference_kinetic_energy=36706.8206525e9,
                                        emittance_norm=(5e-5, 5e-5),
                                        relative_energy_cut=0.001,
                                        seed=1993,
                                        batchMode=True)

    ref_mass = g4man.g4link.getReferenceMass() * 1e9 # convert to eV

    coll1 = g4man.make_xtg4_collimator("TCP.A.B1") # Use new convenience method
    coll2 = g4man.make_xtg4_collimator("TCP.B.B1")

    ## Generate a simple sequence
    # The 1 m drift are split in 2 so the global aperture check at 1 m 
    # can remove low-energy particles with large angles
    line = xt.Line(
        elements=[xt.Drift(length=0.5),
                  xt.Drift(length=0.5),
                  coll1,
                  xt.Multipole(knl=[0, 1.], ksl=[0,0]),
                  xt.Drift(length=0.5),
                  xt.Drift(length=0.5),
                  coll2,
                  xt.Multipole(knl=[0, -1.], ksl=[0,0])
                  ])


    ## Chose a context
    context = xo.ContextCpu()         # For CPU
    # context = xo.ContectCupy()      # For CUDA GPUs
    # context = xo.ContectPyopencl()  # For OpenCL GPUs

    ## Transfer lattice on context and compile tracking code
    line.config.global_xy_limit = 0.1
    line.config.XTRACK_USE_EXACT_DRIFTS = True
    line.build_tracker(_context=context)

    ## Build particle object on context
    n_part = 100
    N_max_expected_products = 50000
    n_turns = 3

    tot_num_part = n_part + N_max_expected_products

    particles = xp.Particles(
        _capacity=tot_num_part,
		p0c = 36900e9,
        mass0 = ref_mass,
        x=np.random.uniform(-1e-2, 1e-2, n_part),
        px=np.random.uniform(-1e-5, 1e-5, n_part),
        y=np.random.uniform(-2e-3, 2e-3, n_part),
        py=np.random.uniform(-3e-5, 3e-5, n_part),
        zeta=np.random.uniform(-1e-2, 1e-2, n_part),
        delta=np.random.uniform(-1e-4, 1e-4, n_part),
        pdg_id=np.full(n_part, xp.pdg.get_pdg_id_from_name('Pb208'))
    )

    ## Track (saving turn-by-turn data)
    print(f"{particles._num_active_particles=}")
    print(f"{particles._num_lost_particles=}")
    print("turn=0")

    for turn in range(n_turns):
        line.track(particles, num_turns=1)

        print(f"{particles._num_active_particles=}")
        print(f"{particles._num_lost_particles=}")
        print(f"{turn=}")

        if particles._num_active_particles==0:
            break

    total_particles = particles._num_active_particles + particles._num_lost_particles
    print(f"Surviving particles: {particles._num_active_particles} / {total_particles}")

    # Allow 1 particle tolerance for rounding errors on different installation
    assert np.isclose(sum(particles.state==1), 1143, atol=1) # Active particles

def main():
    test_fodo()

if __name__ == '__main__':
    main()
