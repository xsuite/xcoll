# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import pytest
try:
    import rpyc
except ImportError as e:
    rpyc = None
import xcoll as xc

if xc.geant4.environment.ready:
    old_bdsim = xc.geant4.environment.bdsim_older_than(compare_version='1.7.7.develop')  # No unstable particles returned in older BDSIM
else:
    old_bdsim = True

engine_params = [
    pytest.param("fluka", marks=pytest.mark.fluka),
    pytest.param("geant4", marks=pytest.mark.geant4)
]

all_engine_params = [
    pytest.param("everest", marks=pytest.mark.everest),
    pytest.param("fluka", marks=pytest.mark.fluka),
    pytest.param("geant4", marks=pytest.mark.geant4)
]

def check_skip(engine, check_old_bdsim=False):
    if engine == "fluka":
        if not xc.fluka.environment.ready:
            pytest.skip("FLUKA installation not found")
    elif engine == "geant4":
        if not xc.geant4.environment.ready:
            pytest.skip("BDSIM+Geant4 installation not found")
        if rpyc is None:
            pytest.skip("rpyc not installed")

def check_skip_old_bdsim(engine, check_old_bdsim=True):
    if engine == "geant4" and check_old_bdsim and old_bdsim:
        pytest.skip("Old BDSIM version detected; skipping tests needing new version")
