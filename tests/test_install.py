# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path
import pytest

import xobjects as xo
import xtrack as xt
import xcoll as xc
from xpart.test_helpers import flaky_assertions, retry
from xobjects.test_helpers import for_all_test_contexts

path = Path(__file__).parent / 'data'


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("beam", [1, 2], ids=["B1", "B2"])
def test_gaps(beam, test_context):
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    coll = xc.BlackAbsorber(length=1.738, angle=127.5)
    name = 'tcp.b6l7.b1' if beam == 1 else 'tcp.b6r7.b2'
    xc.install_elements(line, name, coll, need_apertures=True)
    line.build_tracker()
    tw = line.twiss()
    beta_gamma_rel = line.particle_ref._xobject.gamma0[0]*line.particle_ref._xobject.beta0[0]
    coll.assign_optics(name=name, nemitt_x=3.5e-6, nemitt_y=2.5e-6, twiss=tw, beta_gamma_rel=beta_gamma_rel)

    # compare jaw vs gaps, tilts etc, with align='upstream' and align='downstream'