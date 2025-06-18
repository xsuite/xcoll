# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import pytest
import numpy as np

import xtrack as xt
import xpart as xp
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts

path = xc._pkg_root.parent / 'tests' / 'data'
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize('beam', [1, 2])
def test_line_accessor(beam, test_context):
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    colldb = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=beam)
    assert str(line.collimators) == ''
    assert len(line.collimators) == 0
    colldb.install_everest_collimators(verbose=True, line=line)
    assert str(line.collimators) != ''
    if beam == 1:
        assert len(line.collimators) == 55
    else:
        assert len(line.collimators) == 56
    expected_families = {'tcl4', 'tcl5', 'tcl6', 'tct2', 'tdi', 'tcli', 'tcld', 'tcp3',
                         'tcsg3', 'tcla3', 'tct15', 'tcdq', 'tcsp', 'tcp7', 'tcsg7',
                         'tcla7', 'tct8'}
    assert set(line.collimators.families.keys()) == set(expected_families)
    assert len(line.collimators['tcp7']) == 3
    assert len(line.collimators['tcp3']) == 1
    assert len(line.collimators['tcl4']) == 2
    assert len(line.collimators['tcl5']) == 2
    assert len(line.collimators['tcl6']) == 2
    assert len(line.collimators['tct2']) == 2
    assert len(line.collimators['tdi']) == 3
    assert len(line.collimators['tcli']) == 2
    assert len(line.collimators['tcld']) == 1
    assert len(line.collimators['tcp3']) == 1
    assert len(line.collimators['tcsg3']) == 4
    assert len(line.collimators['tcla3']) == 4
    assert len(line.collimators['tct15']) == 4
    assert len(line.collimators['tcdq']) == 3
    assert len(line.collimators['tcsp']) == 1
    assert len(line.collimators['tcp7']) == 3
    if beam == 1:
        assert len(line.collimators['tcsg7']) == 14
    else:
        assert len(line.collimators['tcsg7']) == 15
    assert len(line.collimators['tcla7']) == 5
    assert len(line.collimators['tct8']) == 2
    for plane in ['d', 'c', 'b']:
        tcp = f"tcp.{plane}6{'l' if beam==1 else 'r'}7.b{beam}"
        assert tcp in line.collimators
        assert tcp in line.collimators['tcp7']
    assert isinstance(line.collimators.gap, dict)
    assert np.isclose(line.collimators['tcp7'].length, 0.6)
    assert np.isclose(line.collimators['tcp7'].gap, 5.0)
    line.collimators['tcp7'].gap = 8.5
    assert np.isclose(line.collimators['tcp7'].gap, 8.5)
    for plane in ['d', 'c', 'b']:
        tcp = f"tcp.{plane}6{'l' if beam==1 else 'r'}7.b{beam}"
        assert np.isclose(line[tcp].gap, 8.5)
    line.collimators['tcp7'].gap = {f"tcp.c6{'l' if beam==1 else 'r'}7.b{beam}": 5.0,
                                    f"tcp.d6{'l' if beam==1 else 'r'}7.b{beam}": 5.0}
    assert isinstance(line.collimators['tcp7'].gap, dict)
    assert len(list(line.collimators['tcp7'].gap.keys())) == 3
    tcp = f"tcp.d6{'l' if beam==1 else 'r'}7.b{beam}"
    assert np.isclose(line.collimators.gap[tcp], 5.0)
    assert np.isclose(line.collimators['tcp7'].gap[tcp], 5.0)
    assert np.isclose(line[tcp].gap, 5.0)
    tcp = f"tcp.c6{'l' if beam==1 else 'r'}7.b{beam}"
    assert np.isclose(line.collimators.gap[tcp], 5.0)
    assert np.isclose(line.collimators['tcp7'].gap[tcp], 5.0)
    assert np.isclose(line[tcp].gap, 5.0)
    tcp = f"tcp.b6{'l' if beam==1 else 'r'}7.b{beam}"
    assert np.isclose(line.collimators.gap[tcp], 8.5)
    assert np.isclose(line.collimators['tcp7'].gap[tcp], 8.5)
    assert np.isclose(line[tcp].gap, 8.5)
    prev_coll = None
    for coll in line.collimators:
        assert isinstance(coll, xc.EverestCollimator)
        assert coll._context.__class__ is test_context.__class__
        assert coll is not prev_coll
        prev_coll = coll
    for name, coll in line.collimators.items():
        assert isinstance(name, str)
        assert isinstance(coll, xc.EverestCollimator)
        assert coll.name == name
        assert name in line.element_names
        assert line.get(name) is coll
