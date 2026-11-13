# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import pytest
import numpy as np
from pathlib import Path

import xobjects as xo
import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts


path = Path(__file__).parent / 'data'


@pytest.mark.gpu
@pytest.mark.xcother
@for_all_test_contexts
@pytest.mark.parametrize('length', [0., 0.2, 0.6, 1.2])
@pytest.mark.parametrize('method', ['direct', 'black_absorbers', 'everest_collimators'])
@pytest.mark.parametrize('apertures_mode', ['single', 'list', 'list_of_pairs'],
                         ids=['single_aperture', 'per_collimator_apertures',
                              'per_collimator_upstream_downstream_pairs'])
def test_install_device_placeholder(length, method, apertures_mode, capsys, test_context):
    # Use test_context from the beginning to fully test installation on
    # the context, even though this is slow and not advised for real use.
    line = xt.Line(elements={
        'before': xt.Drift(length=2., _context=test_context),
        'coll': xt.Device(length=length, _context=test_context),
        'between': xt.Drift(length=2., _context=test_context),
        'coll2': xt.Device(length=length, _context=test_context),
        'after': xt.Drift(length=2., _context=test_context),
    })
    aperture1 = xt.LimitEllipse(a=0.02, b=0.03, _context=test_context)
    aperture2 = xt.LimitEllipse(a=0.05, b=0.07, _context=test_context)
    aperture1_down = xt.LimitEllipse(a=0.021, b=0.031, _context=test_context)
    aperture2_down = xt.LimitEllipse(a=0.051, b=0.071, _context=test_context)
    # 'single': one aperture spec is broadcast to all installed collimators
    #           (used for both their upstream and downstream side).
    # 'list': a list with one aperture spec per collimator assigns a distinct
    #         aperture (for both sides) to each of them.
    # 'list_of_pairs': a list with one [upstream, downstream] pair per
    #         collimator assigns distinct upstream/downstream apertures to
    #         each of them.
    if apertures_mode == 'single':
        apertures = aperture1
        expected = {'coll': (aperture1, aperture1), 'coll2': (aperture1, aperture1)}
    elif apertures_mode == 'list':
        apertures = [aperture1, aperture2]
        expected = {'coll': (aperture1, aperture1), 'coll2': (aperture2, aperture2)}
    else:
        apertures = [[aperture1, aperture1_down], [aperture2, aperture2_down]]
        expected = {'coll': (aperture1, aperture1_down), 'coll2': (aperture2, aperture2_down)}
    if method == 'direct':
        line.xcoll.collimators.install(
            ['coll', 'coll2'],
            [xc.BlackAbsorber(length=0.6, _context=test_context), xc.BlackAbsorber(length=0.6, _context=test_context)],
            apertures=apertures, need_apertures=True)
    else:
        colldb = xc.CollimatorDatabase(
            collimator_dict={'coll': dict(length=0.6, gap=6., material='CFC'),
                             'coll2': dict(length=0.6, gap=6., material='CFC')},
            nemitt_x=3.5e-6, nemitt_y=3.5e-6)
        if method == 'everest_collimators' and isinstance(test_context, xo.ContextPyopencl):
            pytest.skip("Everest not supported on PyOpenCL context.")
        getattr(colldb, f'install_{method}')(line, apertures=apertures)

    expected_class = (xc.EverestCollimator if method == 'everest_collimators'
                      else xc.BlackAbsorber)
    tt = line.get_table()
    for name in ['coll', 'coll2']:
        assert isinstance(line[name], expected_class)
        assert line[name].length == pytest.approx(0.6)
        assert tt['s', f'{name}/aper_start'] == pytest.approx(tt['s_start', name])
        assert tt['s', f'{name}/aper_end'] == pytest.approx(tt['s_end', name])
    assert tt['s_center', 'coll'] == pytest.approx(2. + length / 2)
    assert line.get_length() == pytest.approx(6. + 2 * length)
    assert 'Removed active element' not in capsys.readouterr().out

    # Verify each collimator got the expected upstream/downstream apertures.
    for name in ['coll', 'coll2']:
        exp_up, exp_down = expected[name]
        installed_up = line[f'{name}/aper_start']
        installed_down = line[f'{name}/aper_end']
        assert np.isclose(installed_up.a, exp_up.a)
        assert np.isclose(installed_up.b, exp_up.b)
        assert np.isclose(installed_down.a, exp_down.a)
        assert np.isclose(installed_down.b, exp_down.b)


@pytest.mark.xcother
@pytest.mark.parametrize('method', ['black_absorbers', 'everest_collimators',
                                    'fluka_collimators', 'geant4_collimators'])
def test_install_from_colldb_at_s_center(method):
    # The collimators do not exist in the line yet: their position has to
    # be taken from `s_center` in the CollimatorDatabase.
    line = xt.Line(elements={
        'drift_1': xt.Drift(length=10.),
        'drift_2': xt.Drift(length=10.),
    })
    colldb = xc.CollimatorDatabase(
        collimator_dict={'coll': dict(length=0.6, gap=6., material='CFC', s_center=3.),
                         'coll2': dict(length=0.4, gap=6., material='CFC', s_center=15.)},
        nemitt_x=3.5e-6, nemitt_y=3.5e-6)
    aperture = xt.LimitEllipse(a=0.02, b=0.03)
    getattr(colldb, f'install_{method}')(line, apertures=aperture)

    tt = line.get_table()
    assert tt['s_center', 'coll'] == pytest.approx(3.)
    assert tt['s_center', 'coll2'] == pytest.approx(15.)
    assert line.get_length() == pytest.approx(20.)
    for name in ['coll', 'coll2']:
        assert tt['s', f'{name}/aper_start'] == pytest.approx(tt['s_start', name])
        assert tt['s', f'{name}/aper_end'] == pytest.approx(tt['s_end', name])
        assert line[f'{name}/aper_start'].a == pytest.approx(0.02)
        assert line[f'{name}/aper_end'].b == pytest.approx(0.03)


@pytest.mark.gpu
@pytest.mark.xcother
@for_all_test_contexts
@pytest.mark.parametrize('mode', [None, 'thin', 'thick'])
def test_install_overlapping_devices(mode, capsys, test_context):
    # Use test_context from the beginning to fully test installation on
    # the context, even though this is slow and not advised for real use.
    line = xt.Line(elements={
        'before': xt.Drift(length=2., _context=test_context),
        'instrument': xt.Device(length=0.25, _context=test_context),
        'monitor': xt.Device(length=0.75, _context=test_context),
        'after': xt.Drift(length=2., _context=test_context),
    })
    if mode is not None:
        line.slice_thick_elements([xt.Strategy(
            xt.Uniform(2, mode=mode), element_type=xt.Device)])

    line.xcoll.collimators.install('coll', xc.BlackAbsorber(length=1., _context=test_context), at=2.)

    tt = line.get_table()
    assert tt['s_start', 'coll'] == pytest.approx(2.)
    assert tt['s_end', 'coll'] == pytest.approx(3.)
    assert tt['s_start', 'after'] == pytest.approx(3.)
    assert line.get_length() == pytest.approx(5.)
    assert not any(isinstance(ee, (xt.Device, xt.ThickSliceDevice))
                   for ee in line.elements)
    assert 'Removed active element' not in capsys.readouterr().out


@pytest.mark.xcother
@pytest.mark.parametrize("aper", [None, "auto", "single", "both",
                                  "single_ref_aper", "both_ref_aper"],
                         ids=["without_aper", "auto_aper", "single_aper",
                              "both_aper", "single_ref_aper", "both_ref_aper"])
@pytest.mark.parametrize("beam", [1, 2], ids=["B1", "B2"])
def test_install_single_existing_marker(beam, aper):
    # Installation is a structural, host-side operation. GPU-context
    # installation is covered on small lines above; doing it repeatedly on a
    # 100k-element LHC line causes device-buffer reconstruction to dominate.
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
            aperture = [xt.LimitEllipse(a=0.01, b=0.01),
                        xt.LimitEllipse(a=0.02, b=0.02)]
    line = env[f'lhcb{beam}']
    machine_length = line.get_length()

    # Test absorber
    name = 'tcp.b6l7.b1' if beam == 1 else 'tcp.b6r7.b2'
    assert not isinstance(line[name], xc.BlackAbsorber)
    tt = line.get_table()
    pos_centre = tt['s', name] + line[name].length/2
    coll = xc.BlackAbsorber(length=0.6, angle=127.5)
    line.xcoll.collimators.install(
                            name,
                            coll,
                            apertures=aperture,
                            need_apertures=need_apertures
                        )
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
    coll = xc.EverestCollimator(length=0.6, angle=90,
                                material=xc.materials.MolybdenumGraphite)
    line.xcoll.collimators.install(
                            name,
                            coll,
                            apertures=aperture,
                            need_apertures=need_apertures
                        )
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
@pytest.mark.parametrize("beam", [1, 2], ids=["B1", "B2"])
def test_install_single_no_marker(beam):
    env = xt.load(path / f'sequence_lhc_run3_b{beam}.json')
    line = env[f'lhcb{beam}']
    machine_length = line.get_length()

    # Test absorber
    name = 'test_absorber'
    assert name not in line.element_names
    coll = xc.BlackAbsorber(length=1.738, angle=127.5)
    line.xcoll.collimators.install(
                                name,
                                coll,
                                at=12.4,
                                need_apertures=True,
                                apertures=xt.LimitEllipse(0.4, 0.4)
                            )
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
    el = xc.EverestBlock(length=0.63, material=xc.materials.Silicon)
    line.xcoll.collimators.install(
                                name,
                                el,
                                need_apertures=False,
                                at=17.89
                            )
    assert name in line.element_names
    assert np.isclose(line[name].length, 0.63)
    assert np.isclose(line.get_table()['s', name], 17.89)
    assert isinstance(line[name], xc.EverestBlock)

    # Verify line length did not corrupt
    assert np.isclose(machine_length, line.get_length())


@pytest.mark.xcother
def test_install_reuses_existing_apertures():
    aper_up = xt.LimitEllipse(a=0.02, b=0.03)
    aper_down = xt.LimitEllipse(a=0.04, b=0.05)

    line = xt.Line(elements={
        'before': xt.Drift(length=2.),
        'aper_up': aper_up,
        'coll': xt.Device(length=0.6),
        'aper_down': aper_down,
        'after': xt.Drift(length=2.),
    })

    names_before = list(line.element_names)
    length_before = line.get_length()

    line.xcoll.collimators.install(
        'coll',
        xc.BlackAbsorber(length=0.6),
        need_apertures=True,
    )

    assert line.element_names == names_before
    assert line.get('aper_up') is aper_up
    assert line.get('aper_down') is aper_down
    assert isinstance(line.get('coll'), xc.BlackAbsorber)
    assert np.isclose(line.get_length(), length_before)

    assert 'coll/aper_start' not in line.element_names
    assert 'coll/aper_end' not in line.element_names


@pytest.mark.xcother
@pytest.mark.parametrize("existing_side", ["upstream", "downstream"])
def test_install_reuses_one_existing_aperture(existing_side):
    aperture = xt.LimitEllipse(a=0.02, b=0.03)

    elements = {"before": xt.Drift(length=2.)}
    if existing_side == "upstream":
        elements["existing_aper"] = aperture
    elements["coll"] = xt.Device(length=0.6)
    if existing_side == "downstream":
        elements["existing_aper"] = aperture
    elements["after"] = xt.Drift(length=2.)

    line = xt.Line(elements=elements)

    original_aperture = line.get("existing_aper")
    original_length = line.get_length()
    original_num_elements = len(line.element_names)

    line.xcoll.collimators.install(
        "coll",
        xc.BlackAbsorber(length=0.6),
        need_apertures=True,
    )

    # The collimator was replaced in place.
    assert isinstance(line.get("coll"), xc.BlackAbsorber)
    assert np.isclose(line.get_length(), original_length)

    # Only one additional aperture was installed.
    assert len(line.element_names) == original_num_elements + 1
    assert line.get("existing_aper") is original_aperture

    idx = line.element_names.index("coll")

    if existing_side == "upstream":
        assert line.element_names[idx - 1] == "existing_aper"
        assert line.element_names[idx + 1] == "coll/aper_end"
        assert "coll/aper_start" not in line.element_names
    else:
        assert line.element_names[idx - 1] == "coll/aper_start"
        assert line.element_names[idx + 1] == "existing_aper"
        assert "coll/aper_end" not in line.element_names

    # Both apertures are at the correct longitudinal positions.
    tt = line.get_table()
    assert np.isclose(tt.s[idx - 1], tt.s_start[idx])
    assert np.isclose(tt.s[idx + 1], tt.s_end[idx])

    # The newly installed aperture has the expected geometry.
    new_name = (
        "coll/aper_end" if existing_side == "upstream"
        else "coll/aper_start"
    )
    assert np.isclose(line.get(new_name).a, aperture.a)
    assert np.isclose(line.get(new_name).b, aperture.b)


@pytest.mark.xcother
def test_install_reuses_removed_aperture():
    aperture = xt.LimitEllipse(a=0.02, b=0.03)

    line = xt.Line(elements={
        'before': xt.Drift(length=2.),
        'coll/aper_start': aperture,
        'coll': xt.Marker(),
        'after': xt.Drift(length=2.),
    })

    length_before = line.get_length()
    centre_before = line.get_table()['s', 'coll']

    # The aperture and marker initially occupy the same position.
    assert np.isclose(
        line.get_table()['s', 'coll/aper_start'], centre_before
    )

    # Installing a thick collimator requires the fallback path.
    # The original aperture is removed from the line but remains
    # in the environment.
    line.xcoll.collimators.install(
        'coll',
        xc.BlackAbsorber(length=0.6),
        need_apertures=True,
    )

    assert isinstance(line.get('coll'), xc.BlackAbsorber)
    assert np.isclose(line.get_length(), length_before)

    tt = line.get_table()

    # The collimator remains centred on the original marker.
    assert np.isclose(tt['s_center', 'coll'], centre_before)

    # The original aperture object and name are reused.
    assert line.get('coll/aper_start') is aperture
    assert line.element_names.count('coll/aper_start') == 1

    # Both apertures are installed at the collimator boundaries.
    idx = line.element_names.index('coll')
    assert line.element_names[idx - 1] == 'coll/aper_start'
    assert line.element_names[idx + 1] == 'coll/aper_end'

    assert np.isclose(
        tt['s', 'coll/aper_start'], tt['s_start', 'coll']
    )
    assert np.isclose(
        tt['s', 'coll/aper_end'], tt['s_end', 'coll']
    )

    # The new downstream aperture inherits the original geometry.
    downstream = line.get('coll/aper_end')
    assert np.isclose(downstream.a, aperture.a)
    assert np.isclose(downstream.b, aperture.b)
