import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp

import collimasim as cs

import multiprocessing

def test_fodo():
    np.random.seed(seed=1994)

    g4man = cs.Geant4CollimationManager(#collimator_file="resources/CollDB_new_example.dat",
                                        #collimator_file="resources/CollDB_old_example.dat",
                                        collimator_file="resources/collgaps.dat",
                                        bdsim_config_file="resources/settings.gmad",
                                        tfs_file="resources/collonly_twiss_file_example.tfs",
                                        reference_pdg_id=-11,
                                        reference_kinetic_energy=182.4994890018054e9, # eV
                                        emittance_norm=(0.000521429683675842, 1.0357164949725628e-06),
                                        relative_energy_cut=0.001,
                                        seed=1993,
                                        batchMode=True)

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
    N_max_expected_products = 5000
    n_turns = 3

    tot_num_part = n_part + N_max_expected_products

    particles = xp.Particles(
        _capacity=tot_num_part,
		p0c = 185.e9,
        mass0 = xp.constants.ELECTRON_MASS_EV,
        x=np.random.uniform(-1e-2, 1e-2, n_part),
        px=np.random.uniform(-1e-5, 1e-5, n_part),
        y=np.random.uniform(-2e-3, 2e-3, n_part),
        py=np.random.uniform(-3e-5, 3e-5, n_part),
        zeta=np.random.uniform(-1e-2, 1e-2, n_part),
        delta=np.random.uniform(-1e-4, 1e-4, n_part),
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

    # Allow 2 particle tolerance for rounding errors on different installation
    assert np.isclose(sum(particles.state==1), 8, atol=2) # Active particles

def main():
    test_fodo()

if __name__ == '__main__':
    main()

# The BDSIM link instances are no re-entry safe so can only run 1 test per file
# If multiple instances are required, multiprocessing can offer a solution,
# as spawned processes correctly handle the memory
# Comment out the code below to try this functionality
#  
# def test_object():
#     N_part = 10
#     N_max_expected_products = 500
#     tot_num_part = N_part + N_max_expected_products
#     particles = xp.Particles(
#         _capacity=tot_num_part,
# 		p0c = 185.e9,
#         x = np.zeros(N_part),
# 		y = np.zeros(N_part),
# 		zeta = np.zeros(N_part),
# 		px = np.zeros(N_part),
# 		py = np.zeros(N_part),
# 		delta = np.zeros(N_part))
#
#     particles.x[1] = 0.3
#     particles.x[4] = 0.3
#
#     print(particles.x)
#
#     g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_old_example.dat",
#                                         bdsim_config_file="resources/settings.gmad",
#                                         tfs_file="resources/collonly_twiss_file_example.tfs",
#                                         reference_pdg_id=-11,
#                                         reference_kinetic_energy=182.4994890018054e9, # eV
#                                         emittance_norm=(0.000521429683675842, 1.0357164949725628e-06),
#                                         relative_energy_cut=0.001,
#                                         seed=1993,
#                                         batchMode=True)
#
#     coll1 = xt.BeamInteraction(interaction_process=cs.Geant4Collimator(name="tcp.a.b1",
#                                                                        g4manager=g4man))
#     coll2 = xt.BeamInteraction(interaction_process=cs.Geant4Collimator(name="tcp.a.b1",
#                                                                        g4manager=g4man))
#
#     coll1.track(particles)
#     coll2.track(particles)
#
#     print(f"Surviving after the second collimator: {sum(particles.state == 1)} / {particles.state}")
#     assert sum(particles.state == 1) == 44 # Active particles
#
#
# def main():
#     # Initialising and de-allocating the collimasim extensions
#     # may not clear up the memory properly. This seems to be a Python
#     # limitation. Multiprocessing processes seem to solve this issue 
#     p1 = multiprocessing.Process(target=test_object)
#     p2 = multiprocessing.Process(target=test_fodo)
#
#     p1.start()
#     p1.join()
#
#     p2.start()
#     p2.join()
#
# if __name__ == '__main__':
#     main()
