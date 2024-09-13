# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path
import pytest
import numpy as np

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
def test_install_single_existing_marker(beam, test_context):
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    machine_length = line.get_length()

    # Test absorber
    name = 'tcp.b6l7.b1' if beam == 1 else 'tcp.b6r7.b2'
    assert not isinstance(line[name], xc.BlackAbsorber)
    pos_centre = line.get_s_position(name) + line[name].length/2
    coll = xc.BlackAbsorber(length=0.6, angle=127.5)
    xc.install_elements(line, name, coll, need_apertures=True)
    assert np.isclose(line[name].length, 0.6)
    assert np.isclose(pos_centre - line[name].length/2, line.get_s_position(name))
    assert isinstance(line[name], xc.BlackAbsorber)
    tab = line.get_table()
    idx = tab.mask[[name]][0]
    assert xt.line._is_aperture(line[idx-1], line)
    assert xt.line._is_aperture(line[idx+1], line)

    # Test normal collimator
    name = 'tcp.d6l7.b1' if beam == 1 else 'tcp.d6r7.b2'
    assert not isinstance(line[name], xc.EverestCollimator)
    # We will give the Drift at the location of the collimator a length
    # (and subtract that from the drifts before and after), to test the correct placement
    tab = line.get_table()
    existing_length = 0.12
    line[name].length += existing_length
    idx = tab.mask[[name]][0]
    while True:
        idx -= 1
        if tab.element_type[idx].startswith('Drift'):
            assert line[idx].length > existing_length/2
            line[idx].length -= existing_length/2
            break
    idx = tab.mask[[name]][0]
    while True:
        idx += 1
        if tab.element_type[idx].startswith('Drift'):
            assert line[idx].length > existing_length/2
            line[idx].length -= existing_length/2
            break
    pos_centre = line.get_s_position(name) + line[name].length/2
    coll = xc.EverestCollimator(length=0.6, angle=90, material=xc.materials.MolybdenumGraphite)
    xc.install_elements(line, name, coll, need_apertures=True)
    assert np.isclose(line[name].length, 0.6)
    assert np.isclose(pos_centre - line[name].length/2, line.get_s_position(name))
    assert isinstance(line[name], xc.EverestCollimator)
    tab = line.get_table()
    idx = tab.mask[[name]][0]
    assert xt.line._is_aperture(line[idx-1], line)
    assert xt.line._is_aperture(line[idx+1], line)

    # Test crystal collimator
    name = 'tcpcv.a6l7.b1' if beam == 1 else 'tcpcv.a6r7.b2'
    assert not isinstance(line[name], xc.EverestCrystal)
    pos_centre = line.get_s_position(name) + line[name].length/2
    coll = xc.EverestCrystal(length=0.004, angle=90, lattice='strip', bending_radius=85.10,
                             width=5.0e-3, height=30.0e-3, side='left',
                             material=xc.materials.SiliconCrystal)
    xc.install_elements(line, name, coll, need_apertures=True)
    assert np.isclose(line[name].length, 0.004)
    assert np.isclose(pos_centre - line[name].length/2, line.get_s_position(name))
    assert isinstance(line[name], xc.EverestCrystal)
    tab = line.get_table()
    idx = tab.mask[[name]][0]
    assert xt.line._is_aperture(line[idx-1], line)
    assert xt.line._is_aperture(line[idx+1], line)

    # Test block
    name = 'ip4'
    assert not isinstance(line[name], xc.EverestBlock)
    pos_centre = line.get_s_position(name)
    el = xc.EverestBlock(length=0.63, material=xc.materials.Silicon)
    xc.install_elements(line, name, el, need_apertures=False)
    assert np.isclose(line[name].length, 0.63)
    assert np.isclose(pos_centre - line[name].length/2, line.get_s_position(name))
    assert isinstance(line[name], xc.EverestBlock)

    # Verify line length did not corrupt
    assert np.isclose(machine_length, line.get_length())


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("beam", [1, 2], ids=["B1", "B2"])
def test_install_single_no_marker(beam, test_context):
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    machine_length = line.get_length()

    # Test absorber
    name = 'test_absorber'
    assert name not in line.element_names
    coll = xc.BlackAbsorber(length=1.738, angle=127.5)
    xc.install_elements(line, name, coll, at_s=12.4, need_apertures=True,
                        apertures=xt.LimitEllipse(0.4, 0.4))
    assert name in line.element_names
    assert np.isclose(line[name].length, 1.738)
    assert np.isclose(line.get_s_position(name), 12.4)
    assert isinstance(line[name], xc.BlackAbsorber)
    tab = line.get_table()
    idx = tab.mask[[name]][0]
    assert xt.line._is_aperture(line[idx-1], line)
    assert isinstance(line[idx-1], xt.LimitEllipse)
    assert np.isclose(line[idx-1].a_squ, 0.16)
    assert np.isclose(line[idx-1].b_squ, 0.16)
    assert xt.line._is_aperture(line[idx+1], line)
    assert isinstance(line[idx+1], xt.LimitEllipse)
    assert np.isclose(line[idx+1].a_squ, 0.16)
    assert np.isclose(line[idx+1].b_squ, 0.16)

    # Test block
    name = 'test_block'
    assert name not in line.element_names
    el = xc.EverestBlock(length=0.63, material=xc.materials.Silicon)
    xc.install_elements(line, name, el, need_apertures=False, at_s=17.89)
    assert name in line.element_names
    assert np.isclose(line[name].length, 0.63)
    assert np.isclose(line.get_s_position(name), 17.89)
    assert isinstance(line[name], xc.EverestBlock)

    # Verify line length did not corrupt
    assert np.isclose(machine_length, line.get_length())
