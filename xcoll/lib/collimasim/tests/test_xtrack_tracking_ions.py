import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp

import collimasim as cs

import multiprocessing

def make_particles():
    np.random.seed(seed=1994)
    n_part = 5000
    particles = xp.Particles(
        _capacity=n_part,
		p0c = 36900e9,
        mass0 = 193687676197.13638,
        x=np.random.uniform(-1e-3, 1e-3, n_part),
        px=np.random.uniform(-1e-4, 1e-4, n_part),
        y=np.random.uniform(-1e-3, 1e-3, n_part),
        py=np.random.uniform(-1e-4, 1e-4, n_part),
        zeta=np.random.uniform(-0.1, 0.1, n_part),
        delta=np.random.uniform(-0.4, 0.4, n_part),
    )

    return particles


def run_g4_logitudinal(particles):
    np.random.seed(seed=1994)

    g4man = cs.Geant4CollimationManager(collimator_file="resources/CollDB_testing.dat",
                                        bdsim_config_file="resources/settings_black_absorber.gmad",
                                        tfs_file="resources/twiss_file_testing.tfs",
                                        reference_pdg_id=1000822080,
                                        reference_kinetic_energy=36706820652500.0,
                                        emittance_norm=(5.e-5, 5.e-5),
                                        relative_energy_cut=0.001,
                                        seed=1993,
                                        # batchMode=False
                                        batchMode=True
                                        )

    coll1 = g4man.make_xtg4_collimator("coll_open") # Use new convenience method

    # Generate a simple sequence
    line = xt.Line(
        elements=[coll1,
                  ])

    context = xo.ContextCpu() 
    line.build_tracker(_context=context)

    part_copy = particles.copy()
    line.track(part_copy, num_turns=1)
    part_copy.remove_unused_space()

    return part_copy

def run_drift_logitudinal(particles):

    line = xt.Line(
        elements=[xt.Drift(length=10),
                  ])

    context = xo.ContextCpu() 
    line.config.XTRACK_USE_EXACT_DRIFTS = True
    line.build_tracker(_context=context)

    part_copy = particles.copy()
    line.track(part_copy, num_turns=1)
    part_copy.remove_unused_space()

    return part_copy


def test_xtrack_tracking():
    particles = make_particles()
    part_g4 = run_g4_logitudinal(particles)
    part_dr = run_drift_logitudinal(particles)

    # BDSIM tracks trough some geometry padding now, which is not corrected for now
    # so use a fairly lenient tolerance of 1e-8 for all coordinates over 10 m of drift
    all_close = True
    for var in ['x', 'px', 'y', 'py', 'zeta', 'delta']:
        atol = 1e-6 # For ions, the zeta variable doesn't seem to agree so well, e.g it is up to 1e-7 m off
        rtol = 1e-5
        var_close = np.all(np.isclose(getattr(part_g4, var), getattr(part_dr, var), atol=atol, rtol=rtol))

        all_close &= var_close
        print("Var {}: {}".format(var, var_close))
    
    assert all_close


if __name__ == '__main__':
    test_xtrack_tracking()
