import copy
from at import elements, lattice_pass
from at import atpass, elements
from at.lattice import Lattice
from at import Element, lattice_pass
from at import set_shift, set_tilt

# It is assumed that the pyCollimatorPass module is in pyAT
from at.integrators import pyCollimatorPass

import numpy

def make_particles():
    numpy.random.seed(seed=1994)
    n_part = 1000
    rin = numpy.array([numpy.random.uniform(-1e-3, 1e-3, n_part), 
                       numpy.random.uniform(-1e-5, 1e-5, n_part),
                       numpy.random.uniform(-1e-3, 1e-3, n_part),
                       numpy.random.uniform(-1e-5, 1e-5, n_part),
                       numpy.random.uniform(-1e-2, 1e-2, n_part),
                       #numpy.zeros(n_part),
                       numpy.random.uniform(-1e-4, 1e-4, n_part),
                       #numpy.zeros(n_part),
                       ])

    return rin


def run_particles_g4(rin):
    numpy.random.seed(seed=1994)

    # Initialise the Geant4 interface
    pyCollimatorPass.initialise(bdsim_config_file="resources/settings_black_absorber.gmad",
                                collimator_info="resources/collgaps_pyat_test.dat",
                                reference_pdg_id=-11,
                                reference_Ek=182.4994890018054,
                                relative_energy_cut=0.3,
                                seed=1991)


    # Prepare the collimators in the lattice. The pass method must be specified as
    # pyCollimatorPass and the name must match a collimator defined in BDSIM
    lattice = [
               elements.Element('DUMMY1', Length=10, PassMethod='pyCollimatorPass'),
    ]

    r_original = numpy.copy(rin)
    r_out = lattice_pass(lattice, rin, nturns=1, refpts=range(len(lattice)))

    columns=["x", "px", "y", "py", "delta", "ct"]
    final_values = {}
    for i, var in enumerate(columns):
        #final_values[var] = r_out[i, :, -1, -1]
        final_values[var] = rin[i]

    return final_values

def run_particles_drift(rin):
    lattice = [
               elements.Drift('DRIFT1', 10),
    ]

    r_original = numpy.copy(rin)
    r_out = lattice_pass(lattice, rin, nturns=1, refpts=range(len(lattice)))

    #r_out_dr = lattice_pass(lattice2, r_original, nturns=1, refpts=range(len(lattice2)))

    columns=["x", "px", "y", "py", "delta", "ct"]
    final_values = {}
    for i, var in enumerate(columns):
        #final_values[var] = r_out[i, :, -1, -1]
        final_values[var] = rin[i]

    return final_values

def test_pyat_tracking():
    rin = make_particles()
    rin1 = copy.deepcopy(rin)
    rin2 = copy.deepcopy(rin)

    final_values_g4 = run_particles_g4(rin1)
    final_values_dr = run_particles_drift(rin2)

    # BDSIM tracks trough some geometry padding now, which is not corrected for now
    # so use a bad tolerance of 1e-6 for all coordinates over 10 m of drift
    # TODO: investigate

    all_close = True
    columns=["x", "px", "y", "py", "delta", "ct"]
    for var in columns:
        atol = 1e-6
        rtol = 1e-5
        var_close = numpy.all(numpy.isclose(final_values_dr[var], final_values_g4[var], atol=atol, rtol=rtol))
        all_close &= var_close

        print("Var {}: {}".format(var, var_close))

    assert all_close


if __name__ == "__main__":
    test_pyat_tracking()
