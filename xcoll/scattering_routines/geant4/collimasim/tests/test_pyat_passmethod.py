from at import elements, lattice_pass
from at import atpass, elements
from at.lattice import Lattice
from at import Element, lattice_pass
from at import set_shift, set_tilt

# It is assumed that the pyCollimatorPass module is in pyAT
from at.integrators import pyCollimatorPass

import numpy


def test_multiple_particles_collimator_pass():
    # Initialise the Geant4 interface
    pyCollimatorPass.initialise(bdsim_config_file="resources/settings.gmad",
                                collimator_info="resources/collgaps.dat",
                                reference_pdg_id=-11,
                                reference_Ek=182.4994890018054,
                                relative_energy_cut=0.3,
                                seed=1991)


    # Prepare the collimators in the lattice. The pass method must be specified as
    # pyCollimatorPass and the name must match a collimator defined in BDSIM
    lattice = [
               elements.Drift('drift_upstream', 0.6, PassMethod='pyDriftPass'),
               elements.Element('TCP.A.B1', Length=0.6, PassMethod='pyCollimatorPass'),
               elements.Element('TCP.B.B1', Length=0.6, PassMethod='pyCollimatorPass'),
               elements.Drift('drift_downstream', 0.6)
    ]

    # Track 5 test particles
    rin = numpy.zeros((6, 5))

    rin[0, 0] = 9.762728e-3  # particle one offset in x
    rin[2, 1] = 10e-3  # particle two offset in y
    rin[2, 2] = 2e-3
    rin[2, 3] = 0.5e-3
    # The 5th particle is on the reference orbit

    r_original = numpy.copy(rin)
    r_out = lattice_pass(lattice, rin, nturns=1, refpts=range(len(lattice)))

    columns=["x", "px", "y", "py", "delta", "ct"]
    final_values = {}
    for i, var in enumerate(columns):
        final_values[var] = r_out[i, :, -1, -1]

    print("final values:\n")
    per_particle = list(zip(final_values['x'], final_values['y'], final_values['px'], final_values['py'], final_values['delta'], final_values['ct']))
    for pp in per_particle:
        print(pp)

    assert(sum(numpy.isnan(r_out[0, :, -1, -1])) == 2) # Number of lost partciles


if __name__ == "__main__":
    test_multiple_particles_collimator_pass()

