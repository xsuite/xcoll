# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import pytest
import numpy as np
from pathlib import Path

import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts


path = Path(__file__).parent / 'data'


@pytest.mark.xcother
@pytest.mark.parametrize('length', [0., 0.2, 1.2])
@pytest.mark.parametrize('method', ['direct', 'black_absorbers', 'everest_collimators'])
def test_install_device_placeholder(length, method, capsys):
    line = xt.Line(elements={
        'before': xt.Drift(length=2.),
        'coll': xt.Device(length=length),
        'after': xt.Drift(length=2.),
    })
    aperture = xt.LimitEllipse(a=0.02, b=0.03)
    if method == 'direct':
        line.xcoll.collimators.install(
            'coll', xc.BlackAbsorber(length=0.6),
            apertures=aperture, need_apertures=True)
    else:
        colldb = xc.CollimatorDatabase(
            collimator_dict={'coll': dict(length=0.6, gap=6., material='CFC')},
            nemitt_x=3.5e-6, nemitt_y=3.5e-6)
        getattr(colldb, f'install_{method}')(line, apertures=aperture)

    expected_class = (xc.EverestCollimator if method == 'everest_collimators'
                      else xc.BlackAbsorber)
    assert isinstance(line['coll'], expected_class)
    assert line['coll'].length == pytest.approx(0.6)
    tt = line.get_table()
    assert tt['s_center', 'coll'] == pytest.approx(2. + length / 2)
    assert tt['s', 'coll_aper_upstream'] == pytest.approx(tt['s_start', 'coll'])
    assert tt['s', 'coll_aper_downstream'] == pytest.approx(tt['s_end', 'coll'])
    assert line.get_length() == pytest.approx(4. + length)
    assert 'Removed active element' not in capsys.readouterr().out


@pytest.mark.xcother
@pytest.mark.parametrize('mode', [None, 'thin', 'thick'])
def test_install_overlapping_devices(mode, capsys):
    line = xt.Line(elements={
        'before': xt.Drift(length=2.),
        'instrument': xt.Device(length=0.25),
        'monitor': xt.Device(length=0.75),
        'after': xt.Drift(length=2.),
    })
    if mode is not None:
        line.slice_thick_elements([xt.Strategy(
            xt.Uniform(2, mode=mode), element_type=xt.Device)])

    line.xcoll.collimators.install('coll', xc.BlackAbsorber(length=1.), at=2.)

    tt = line.get_table()
    assert tt['s_start', 'coll'] == pytest.approx(2.)
    assert tt['s_end', 'coll'] == pytest.approx(3.)
    assert tt['s_start', 'after'] == pytest.approx(3.)
    assert line.get_length() == pytest.approx(5.)
    assert not any(isinstance(ee, (xt.Device, xt.ThickSliceDevice))
                   for ee in line.elements)
    assert 'Removed active element' not in capsys.readouterr().out


@pytest.mark.xcother
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("aper", [None, "auto", "single", "both",
                                  "single_ref_aper", "both_ref_aper"],
                         ids=["without_aper", "auto_aper", "single_aper",
                              "both_aper", "single_ref_aper", "both_ref_aper"])
@pytest.mark.parametrize("beam", [1, 2], ids=["B1", "B2"])
def test_install_single_existing_marker(beam, aper, test_context):
    aperture = None
    need_apertures = aper is not None
    if aper == 'auto' or (aper is not None and aper.endswith('_ref_aper')):
        env = xt.load(path / f'sequence_lhc_run3_b{beam}.json')
        if aper == 'single_ref_aper':
            aperture = 'tcp.b6l7.b1_aper' if beam == 1 else 'tcp.b6r7.b2_aper'
        elif aper == 'both_ref_aper':
            if beam == 1:
                aperture = ['tcp.b6l7.b1_aper', 'tcp.d6l7.b1_aper']
            else:
                aperture = ['tcp.b6r7.b2_aper', 'tcp.d6r7.b2_aper']
    else:
        env = xt.load(path / f'sequence_lhc_run3_b{beam}_no_aper.json')
        if aper == 'single':
            aperture = xt.LimitEllipse(a=0.01, b=0.01)
        elif aper == 'both':
            aperture = [xt.LimitEllipse(a=0.01, b=0.01), xt.LimitEllipse(a=0.02, b=0.02)]
    line = env[f'lhcb{beam}']
    machine_length = line.get_length()

    # Test absorber
    name = 'tcp.b6l7.b1' if beam == 1 else 'tcp.b6r7.b2'
    assert not isinstance(line[name], xc.BlackAbsorber)
    tt = line.get_table()
    pos_centre = tt['s', name] + line[name].length/2
    coll = xc.BlackAbsorber(length=0.6, angle=127.5, _context=test_context)
    line.xcoll.collimators.install(name, coll, apertures=aperture, need_apertures=need_apertures)
    assert np.isclose(line[name].length, 0.6)
    assert np.isclose(pos_centre - line[name].length/2,
                      line.get_table()['s', name])
    assert isinstance(line[name], xc.BlackAbsorber)
    tab = line.get_table()
    if need_apertures:
        idx = tab.rows.indices[[name]][0]
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
    idx = tab.rows.indices[[name]][0]
    while True:
        idx -= 1
        if tab.element_type[idx].startswith('Drift'):
            assert line[idx].length > existing_length/2
            line[idx].length -= existing_length/2
            break
    idx = tab.rows.indices[[name]][0]
    while True:
        idx += 1
        if tab.element_type[idx].startswith('Drift'):
            assert line[idx].length > existing_length/2
            line[idx].length -= existing_length/2
            break
    tt = line.get_table()
    pos_centre = tt['s', name] + line[name].length/2
    coll = xc.EverestCollimator(length=0.6, angle=90, _context=test_context,
                                material=xc.materials.MolybdenumGraphite)
    line.xcoll.collimators.install(name, coll, apertures=aperture, need_apertures=need_apertures)
    assert np.isclose(line[name].length, 0.6)
    assert np.isclose(pos_centre - line[name].length/2,
                      line.get_table()['s', name])
    assert isinstance(line[name], xc.EverestCollimator)
    tab = line.get_table()
    if need_apertures:
        idx = tab.rows.indices[[name]][0]
        assert xt.line._is_aperture(line[idx-1], line)
        assert xt.line._is_aperture(line[idx+1], line)

    # Verify line length did not corrupt
    assert np.isclose(machine_length, line.get_length())


@pytest.mark.xcother
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("beam", [1, 2], ids=["B1", "B2"])
def test_install_single_no_marker(beam, test_context):
    env = xt.load(path / f'sequence_lhc_run3_b{beam}.json')
    line = env[f'lhcb{beam}']
    machine_length = line.get_length()

    # Test absorber
    name = 'test_absorber'
    assert name not in line.element_names
    coll = xc.BlackAbsorber(length=1.738, angle=127.5, _context=test_context)
    line.xcoll.collimators.install(name, coll, at=12.4, need_apertures=True,
                             apertures=xt.LimitEllipse(0.4, 0.4))
    assert name in line.element_names
    assert np.isclose(line[name].length, 1.738)
    assert np.isclose(line.get_table()['s', name], 12.4)
    assert isinstance(line[name], xc.BlackAbsorber)
    tab = line.get_table()
    idx = tab.rows.indices[[name]][0]
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
    el = xc.EverestBlock(length=0.63, material=xc.materials.Silicon, _context=test_context)
    line.xcoll.collimators.install(name, el, need_apertures=False, at=17.89)
    assert name in line.element_names
    assert np.isclose(line[name].length, 0.63)
    assert np.isclose(line.get_table()['s', name], 17.89)
    assert isinstance(line[name], xc.EverestBlock)

    # Verify line length did not corrupt
    assert np.isclose(machine_length, line.get_length())
