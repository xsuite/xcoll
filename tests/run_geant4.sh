#!/bin/bash
pytest test_geant4.py::test_serial_bdsim
pytest --html="pytest_geant4_results.html" -n 12 test_geant4.py test_lossmap.py::test_lossmap_geant4 test_jaw_position.py::test_geant4 test_particle_id_hierarchy.py test_ions.py
