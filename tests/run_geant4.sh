#!/bin/bash
pytest test_geant4_serial_bdsim.py
pytest --html="pytest_geant4_results.html" -vv -n 12 test_geant4.py test_geant4_many_particles.py test_lossmap.py::test_lossmap_geant4 test_jaw_position.py::test_geant4 test_particle_id_hierarchy.py test_ions.py
