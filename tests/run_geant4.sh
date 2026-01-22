#!/bin/bash
python ../examples/geant4_init.py
if [ $? -ne 0 ]; then
    echo "Geant4 initialization failed"
    exit 1
fi
pytest test_geant4_serial_bdsim.py
if [ $? -ne 0 ]; then
    echo "Geant4 serial test failed"
    exit 1
fi
pytest --html="pytest_geant4_results.html" -vv -n 12 -m geant4
