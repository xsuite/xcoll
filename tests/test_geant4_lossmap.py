# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
import numpy as np
from pathlib import Path
import xpart as xp
import xtrack as xt
import xcoll as xc
import pytest
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts


# try the import here and skip tests if missing
# also need the import here in case of pytest --forked
try:
    import collimasim as cs
except ImportError:
    cs = None   

#path = Path(__file__).parent / 'data'
path = Path('/afs/cern.ch/work/b/bjlindst/public/git/xtrack/dev/xcoll_geant4/xcoll/tests/data')


beam = 1
plane = 'H'
npart = 500
interpolation = 0.1

line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')

colldb = xc.CollimatorDatabase.from_yaml(path / f'colldb_lhc_run3_ir7.yaml',
                                    beam=beam, ignore_crystals=True)

colldb.install_geant4_collimators(verbose=False, line=line, random_seed=1993,
                                      bdsim_config_file=str(path / f'../data_test_geant4/settings_protons.gmad'))

df_with_coll = line.check_aperture()
assert not np.any(df_with_coll.has_aperture_problem)

line.build_tracker()
xc.assign_optics_to_collimators(line=line)

tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"

part = xc.generate_pencil_on_collimator(line, tcp, num_particles=npart, _capacity=3*npart)

xc.enable_scattering(line)
line.track(part, num_turns=2)
xc.disable_scattering(line)

line_is_reversed = True if beam == 2 else False
ThisLM = xc.LossMap(line, line_is_reversed=line_is_reversed, part=part,
                         interpolation=interpolation)
ThisLM.to_json("lossmap.json")

