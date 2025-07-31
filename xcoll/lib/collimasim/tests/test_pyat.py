import numpy as np

import collimasim


def test_pyat_object():
    # Initialise a code interface - in this case the interface to pyAT
    g4link = collimasim.PyATInterface(bdsimConfigFile="resources/settings.gmad", referencePdgId=-11, referenceEk=182.4994890018054,
                                    relativeEnergyCut=0.001, seed=1993)

    # Add collimators to the Geant4 model
    # Collimators are made of block jaws and are prepared individually in cells, surrounded by a perfect absorber
    # When a collimator is selected, the beam particles start at the beginning of the cell and are extracted at the end
    g4link.addCollimator(name="TESTCOLL1", material="C", length=0.6, aperture=0.1, rotation=0, xOffset=0, yOffset=0, jawTiltLeft=0, jawTiltRight=0, side=0)
    g4link.addCollimator("TESTCOLL2", "C", 0.6, 0.1, np.pi/2, 0, 0, 0, 0, 0)


    # Define some test particles
    #r = np.array([0.1, 0.001, 0.1, 0.001, 0, 0])
    r = np.array([0.0, 0., 0.0, 0., 0., 0.])
    r2 = np.array([0.051, 0., 0.051, 0., 0., 0.])
    r3 = np.array([0.0, 0., 0., 0., 0., 0.])

    # Add the particles to the model
    g4link.addParticle(r)
    g4link.addParticle(r2)
    g4link.addParticle(r3)

    # Select a collimator
    g4link.selectCollimator("TESTCOLL1")

    # Perform the Geant4 simulation
    g4link.collimate()

    # Collect the surviving particles - an n x 6 array where n is the number of surviving particles
    arr_out = g4link.collimateReturn()

    # For lost particles, all coordinates are NaN
    # print(arr_out) # Debug
    print(f"Pass {1}, Output size: {int(len(arr_out)/6)}, Active particles: {int(np.count_nonzero(~np.isnan(arr_out))/6)}")

    # Repeat the collimator pass a few times
    for k in range(10):
        g4link.clearData()

        # Add the particles to the model
        for i in range(0, len(arr_out), 6):
            g4link.addParticle(arr_out[i:i+6])

        # Select a collimator
        g4link.selectCollimator("TESTCOLL1")

        # Perform the Geant4 simulation
        g4link.collimate()

        # Collect the surviving particles - an n x 6 array where n is the number of surviving particles
        arr_out = g4link.collimateReturn()

        # print(arr_out) # Debug
        print(f"Pass {k+2}, Output size: {int(len(arr_out)/6)}, Active particles: {int(np.count_nonzero(~np.isnan(arr_out))/6)}")

    assert int(np.count_nonzero(~np.isnan(arr_out))/6) == 1

if __name__ == "__main__":
    test_pyat_object()