# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import warnings

import numpy as np
import pytest
import xobjects as xo
import xtrack as xt
import xcoll as xc
from xcoll import initial_distribution
from xcoll.optics import _get_optics_at_positions


@pytest.fixture
def make_line():
    def make():
        line = xt.Line(elements={
            'before': xt.Drift(length=2.),
            'coll': xt.Device(length=0.6),
            'between': xt.Drift(length=2.),
            'crystal': xt.Device(length=0.004),
            'after': xt.Drift(length=2.),
        }, particle_ref=xt.Particles(p0c=7e12))
        line.xcoll.collimators.install(
            ['coll', 'crystal'],
            [xc.EverestCollimator(length=0.6, material=xc.materials.MolybdenumGraphite,
                                jaw=0.002),
             xc.EverestCrystal(length=0.004, material=xc.materials.Silicon,
                               jaw=0.002, bending_angle=50e-6)],
            apertures=xt.LimitEllipse(a=0.04, b=0.04), need_apertures=True)
        return line
    return make


@pytest.mark.xcother
def test_installed_element_copy(make_line, monkeypatch):
    line = make_line()
    original_copy = xt.BeamElement.copy
    originals = [line.get(name) for name in ('coll', 'crystal')]

    def check_copy(source, **kwargs):
        # Only the element being copied is temporarily detached.
        # Nested BeamElements (e.g. RandomRutherford) are also copied.
        if any(source is el for el in originals):
            assert '_line' not in source.__dict__
            assert all(el._line is line for el in originals if el is not source)
        return original_copy(source, **kwargs)

    def forbid_line_copy(*args, **kwargs):
        pytest.fail('Copying a single element must not copy its line')

    with monkeypatch.context() as m:
        m.setattr(xt.Line, 'copy', forbid_line_copy)
        m.setattr(xt.BeamElement, 'copy', check_copy)
        for raw in originals:
            copied = raw.copy(_context=xo.ContextCpu())
            assert type(copied) is type(raw)
            assert copied is not raw
            assert copied.name == raw.name
            assert copied.length == raw.length
            assert copied.jaw == raw.jaw
            assert copied._context is not raw._context
            assert copied._line is line
            assert raw._line is line
            copied.jaw = 0.003
            assert raw.jaw != copied.jaw

    def fail_copy(source, **kwargs):
        assert source is originals[0]
        assert '_line' not in source.__dict__
        raise RuntimeError('copy failed')

    with monkeypatch.context() as m:
        m.setattr(xt.BeamElement, 'copy', fail_copy)
        with pytest.raises(RuntimeError, match='copy failed'):
            originals[0].copy()
        assert originals[0]._line is line

    # An unattached element still follows the ordinary copy path.
    raw = xc.BlackAbsorber(length=0.6, jaw=0.002)
    assert '_line' not in raw.copy().__dict__


@pytest.mark.xcother
def test_line_copy_uses_its_own_distribution_api(make_line, monkeypatch):
    line = make_line()
    copied = line.copy(_context=xo.ContextCpu())
    shallow = line.copy(shallow=True)
    for name in ('coll', 'crystal'):
        assert copied.get(name) is not line.get(name)
        # Transitional compatibility: the legacy pointer still refers to the
        # original line. The new API must ignore it.
        assert copied.get(name)._line is line
        assert line.get(name)._line is line
        assert shallow.get(name) is line.get(name)

    def pencil(**kwargs):
        return kwargs

    monkeypatch.setattr(initial_distribution, 'generate_pencil_on_collimator', pencil)
    for current in (line, copied, shallow):
        for name in ('coll', 'crystal'):
            # The explicit API must not depend on the element's old line pointer.
            result = current.xcoll.collimators.generate_pencil(name, 17)
            assert result['line'] is current
            assert result['name'] == name
            assert result['num_particles'] == 17


@pytest.mark.xcother
def test_distribution_api_and_deprecations(make_line, monkeypatch):
    line = make_line()
    api = line.xcoll.collimators
    monkeypatch.setattr(initial_distribution, 'generate_pencil_on_collimator',
                        lambda **kwargs: kwargs)
    monkeypatch.setattr(initial_distribution, 'generate_delta_from_dispersion',
                        lambda **kwargs: kwargs)
    twiss = object()
    tw = object()
    options = dict(side='+', pencil_spread=2e-6, impact_parameter=3e-6,
                   sigma_z=0.08, twiss=twiss, longitudinal='bucket',
                   longitudinal_betatron_cut=2, tw=tw, _capacity=32)
    # Test the complete forwarding contract, including the legacy tw alias.
    # The mock intentionally does not evaluate twiss or the random distribution.
    with warnings.catch_warnings():
        warnings.filterwarnings('error', category=FutureWarning)
        for name in ('coll', 'crystal'):
            result = api.generate_pencil(name, 10, **options)
            assert result['line'] is line
            assert result['name'] == name
            assert result['num_particles'] == 10
            assert all(result[key] == val for key, val in options.items())
        delta_options = dict(plane='x', position_mm=0.002, nemitt_x=2.5e-6,
                             nemitt_y=3.5e-6, betatron_cut=1,
                             match_at_front=False, twiss=twiss)
        delta = api.generate_delta('coll', **delta_options)
        assert delta['line'] is line
        assert delta['at_element'] == 'coll'
        assert all(delta[key] == val for key, val in delta_options.items())

    # Only these explicit compatibility tests exercise the deprecated API.
    for name in ('coll', 'crystal'):
        with pytest.warns(FutureWarning, match=r'Element.generate_pencil') as record:
            result = line[name].generate_pencil(num_particles=10, **options)
        assert len(record) == 1
        assert result['line'] is line and result['name'] == name
        assert all(result[key] == val for key, val in options.items())
    with pytest.warns(FutureWarning, match=r'Element.generate_delta'):
        result = line['coll'].generate_delta(**delta_options)
    assert result == delta

    # Preserve the old association for the transitional element API.
    copied = line.get('coll').copy()
    with pytest.warns(FutureWarning, match=r'Element.generate_pencil'):
        result = copied.generate_pencil(10, **options)
    assert result['line'] is line and result['name'] == 'coll'

    # A genuinely standalone element has no line association.
    detached = xc.BlackAbsorber(length=0.6, jaw=0.002)
    with pytest.warns(FutureWarning, match=r'Element.generate_pencil'):
        with pytest.raises(ValueError, match='missing a pointer'):
            detached.generate_pencil(10)


@pytest.mark.xcother
def test_optics_api_and_element_methods(make_line, monkeypatch):
    line = make_line()
    api = line.xcoll.collimators
    # Use a precomputed table to exercise the actual row selection and optics
    # assignment without a tracking kernel or a periodic optics solution.
    table = line.get_table()
    n = len(table.name)
    twiss = xt.Table(dict(name=table.name.copy(), s=table.s.copy(),
                          betx=np.full(n, 10.), bety=np.full(n, 12.),
                          alfx=np.full(n, -0.5), alfy=np.zeros(n),
                          x=np.zeros(n), y=np.zeros(n)))
    twiss_calls = []

    def compute_twiss(self):
        twiss_calls.append(self)
        return twiss

    monkeypatch.setattr(xt.Line, 'twiss', compute_twiss)
    coll = line.get('coll')
    crystal = line.get('crystal')
    with warnings.catch_warnings():
        warnings.filterwarnings('error', category=FutureWarning)
        # A supplied table avoids a new Twiss computation.
        api.assign_optics(names='coll', nemitt_x=2.5e-6, nemitt_y=3.5e-6,
                          twiss=twiss)
        assert coll.optics_ready()
        assert len(crystal.reference_optics) == 0
        assert not crystal.optics_ready()
        assert not twiss_calls
        # With no table, compute Twiss only once for both collimators.
        api.assign_optics(nemitt_x=2.5e-6, nemitt_y=3.5e-6)
        assert len(twiss_calls) == 1 and twiss_calls[0] is line
        for name in ('coll', 'crystal'):
            element = line.get(name)
            assert element.optics_ready()
            assert isinstance(element.reference_optics, xt.Table)
            assert list(element.reference_optics.position) == ['upstream', 'centre', 'downstream']
            np.testing.assert_allclose(element.reference_optics.s_rel, [0., element.length/2, element.length])
            assert element.reference_optics['s_rel', 'centre'] == pytest.approx(element.length/2)
            assert element.reference_optics.collimator_name == name
            assert element.reference_optics.s_start == pytest.approx(twiss['s', name])
            assert element.reference_optics['betx', 'upstream'] == twiss['betx', name]

        # The public element methods remain useful with explicit optics,
        # including when the element is not installed in a line.
        standalone = xc.BlackAbsorber(length=coll.length, jaw=0.002)
        standalone.assign_optics(name='coll', twiss=twiss, beta_gamma_rel=1.,
                                 nemitt_x=2.5e-6, nemitt_y=3.5e-6)
        assert standalone.optics_ready()
        assert '_line' not in standalone.__dict__
        optics = api.get_optics_at('crystal', positions=['upstream', 'centre', 'downstream'],
                                   twiss=twiss)
        standalone_crystal = xc.BlackCrystal(length=crystal.length, gap=6.)
        standalone_crystal.assign_optics(twiss=optics, beta_gamma_rel=1.,
                                         nemitt_x=2.5e-6, nemitt_y=3.5e-6)
        standalone_crystal.align_to_beam_divergence()
        assert standalone_crystal.tilt == pytest.approx(
            standalone_crystal.divergence * standalone_crystal.gap, abs=5.1e-11)

        crystal.gap = 6.
        api.align_to_beam_divergence(names='crystal')
        expected_tilt = crystal.divergence * crystal.gap
        assert crystal.tilt == pytest.approx(expected_tilt, abs=5.1e-11)
        crystal.tilt = 0.
        api.align_to_beam_divergence()
        assert crystal.tilt == pytest.approx(expected_tilt, abs=5.1e-11)
        assert len(twiss_calls) == 1
        # Extra interior positions are supported, and do not change jaw/gap.
        jaw = coll.jaw
        api.assign_optics(names='coll', twiss=twiss,
                          positions=['upstream', 0.15, 'centre', 'downstream'])
        assert coll.reference_optics['s_rel', 's=0.15'] == pytest.approx(0.15)
        assert list(coll.reference_optics.rows[0.1:0.4:'s_rel'].position) == ['s=0.15', 'centre']
        assert coll.jaw == jaw
        assert len(twiss_calls) == 1
    # Backwards-compatible tuple and deprecated element arguments.
    with pytest.warns(FutureWarning, match='optics'):
        old_row = coll.reference_optics['upstream']
    assert old_row.s_rel[0] == 0.
    assert old_row.s[0] == twiss['s', 'coll']
    batch = api.get_optics_at(['coll', 'crystal'], positions=['centre'], twiss=twiss)
    assert isinstance(batch, xt.Table)
    assert set(batch.name) == {'coll', 'crystal'}
    assert batch.rows[batch.name == 'crystal'].s_rel[0] == pytest.approx(crystal.length/2)
    assigned = api.reference_optics
    assert set(assigned.name) == {'coll', 'crystal'}
    assert assigned.nemitt_x == 2.5e-6
    centre = assigned.rows[(assigned.name == 'coll') & (assigned.position == 'centre')]
    assert len(centre) == 1
    assert centre.s_rel[0] == pytest.approx(coll.length/2)
    np.testing.assert_allclose(assigned.rows[assigned.name == 'coll'].s_rel,
                               coll.reference_optics.s_rel)
    assigned.betx[:] = -1.
    assert coll.reference_optics['betx', 'centre'] > 0.
    coll.set_reference_emittance(nemitt_x=4e-6)
    assert coll.reference_optics.nemitt_x == 4e-6
    assigned = api.reference_optics
    assert 'nemitt_x' in assigned._col_names
    assert np.all(assigned.rows[assigned.name == 'coll'].nemitt_x == 4e-6)
    standalone.assign_optics(reference_optics=assigned, name='coll')
    assert standalone.reference_optics.nemitt_x == 4e-6
    api.assign_optics(reference_optics=assigned, names='coll')
    assert len(twiss_calls) == 1
    assert coll.reference_optics.nemitt_x == 4e-6
    local = api.get_optics_at('coll', positions=['downstream', 'upstream', 'centre'], twiss=twiss)
    assert isinstance(local, xt.Table)
    assert list(local.position) == ['upstream', 'centre', 'downstream']
    assert local.s[0] == pytest.approx(twiss['s', 'coll'])
    assert 'beta_gamma_rel' not in local.keys()
    standalone.assign_optics(twiss=local, beta_gamma_rel=2.)
    assert standalone.reference_optics.beta_gamma_rel == 2.
    assert 'beta_gamma_rel' not in local.keys()
    assert standalone.reference_optics is not local
    copied = standalone.copy()
    assert copied.reference_optics is not standalone.reference_optics
    np.testing.assert_allclose(copied.reference_optics.s, standalone.reference_optics.s)
    assert copied.reference_optics.beta_gamma_rel == standalone.reference_optics.beta_gamma_rel
    for name, element in [('coll', standalone), ('crystal', standalone_crystal)]:
        up, down = api.get_optics_at(name, twiss=twiss)
        with pytest.warns(FutureWarning, match='twiss_upstream') as record:
            element.assign_optics(twiss_upstream=up, twiss_downstream=down,
                                  beta_gamma_rel=1., nemitt_x=2.5e-6, nemitt_y=3.5e-6)
        assert len(record) == 1
        assert element.optics_ready()
    with pytest.raises(ValueError, match='valid crystal'):
        api.align_to_beam_divergence(names='coll')
    with pytest.raises(ValueError, match='valid collimator'):
        api.assign_optics(names='before', twiss=twiss)


@pytest.mark.xcother
def test_optics_at_positions_drift():
    # Compare against a real Xtrack drift with a marker at its centre.
    # The two halves are pure drifts, so no Xcoll slicing is involved.
    line = xt.Line(elements={'d1': xt.Drift(length=0.3),
                            'mid': xt.Marker(),
                            'd2': xt.Drift(length=0.3),
                            'end': xt.Marker()},
                   particle_ref=xt.Particles(p0c=7e12))
    tw = line.twiss(method='4d', betx=10., bety=12., alfx=-0.5, alfy=0.2,
                    x=0.001, px=0.0001, y=-0.002, py=-0.0002,
                    dx=0.1, dpx=0.02, dy=-0.03, dpy=0.01)
    # Remove the centre rows to emulate a normal element-boundary Twiss table.
    boundary = tw.rows[['d1', 'end']]
    optics = _get_optics_at_positions(boundary, ['d1'], [0.6],
                                     ['upstream', 'centre', 'downstream', 0.3])
    for key in ('centre', 's=0.3'):
        row = optics.rows[optics.position == key]
        assert row.s[0] == pytest.approx(tw['s', 'mid'])
        for field in ('betx', 'bety', 'alfx', 'alfy', 'gamx', 'gamy',
                      'x', 'y', 'dx', 'dy', 'px', 'py', 'dpx', 'dpy'):
            # Small-angle drift model; Xtrack can use the exact drift map.
            assert row[field][0] == pytest.approx(tw[field, 'mid'], rel=1e-6, abs=1e-10)
        assert 'W_matrix' not in row._col_names
    with pytest.raises(ValueError, match='outside'):
        _get_optics_at_positions(boundary, ['d1'], [0.6], [0.7])
    with pytest.raises(ValueError, match='Invalid optics position'):
        _get_optics_at_positions(boundary, ['d1'], [0.6], ['middle'])
    with pytest.raises(ValueError, match='not compatible'):
        _get_optics_at_positions(boundary, ['d1'], [0.5], ['centre'])


@pytest.mark.xcother
@pytest.mark.parametrize('cls', [xc.BlackAbsorber, xc.BlackCrystal], ids=['collimator', 'crystal'])
def test_reference_emittance(cls):
    with warnings.catch_warnings():
        warnings.filterwarnings('error', category=FutureWarning)
        coll = cls(length=0.6, gap=6., emittance=2.5e-6)
        ref = coll.reference_optics
        assert isinstance(ref, xt.Table)
        assert len(ref) == 0
        assert ref.nemitt_x == ref.nemitt_y == 2.5e-6
        assert ref.gemitt_x is None
        assert not coll.optics_ready()
        assert coll.optics is None
        assert '_nemitt_x' not in coll.__dict__ and '_nemitt_y' not in coll.__dict__
        coll.set_reference_emittance(nemitt_y=3.5e-6)
        assert ref.nemitt_y == 3.5e-6
        assert not coll.optics_ready()
        dct = coll.to_dict()
        restored = cls.from_dict(dct)
        assert restored.reference_optics.nemitt_x == 2.5e-6
        assert restored.reference_optics.nemitt_y == 3.5e-6
        assert len(restored.reference_optics) == 0
        if cls is xc.BlackCrystal:
            # Both bending fields remain in the serialised format. Reading
            # them must preserve the geometry and not modify the input.
            assert 'bending_radius' in dct and 'bending_angle' in dct
            assert restored.bending_radius == pytest.approx(coll.bending_radius)
            assert restored.bending_angle == pytest.approx(coll.bending_angle)
            for radius in (-2., np.inf):
                crystal = cls(length=0.6, gap=6., bending_radius=radius)
                saved = crystal.to_dict()
                clone = cls.from_dict(saved)
                assert clone.bending_radius == pytest.approx(radius)
                assert clone.bending_angle == pytest.approx(crystal.bending_angle)
                assert 'bending_angle' in saved
            bad = dict(dct, bending_angle=0.1)
            with pytest.raises(ValueError, match='Inconsistent bending_radius'):
                cls.from_dict(bad)
            with pytest.raises(ValueError, match='Need to choose between'):
                cls(length=0.6, bending_radius=1., bending_angle=0.1)
        # Both planes are validated before either is changed.
        with pytest.raises(ValueError, match='positive'):
            coll.set_reference_emittance(nemitt_x=4e-6, nemitt_y=-1.)
        assert ref.nemitt_x == 2.5e-6
        with pytest.raises(AttributeError, match='read-only'):
            ref.nemitt_x = 4e-6
        with pytest.raises(AttributeError, match='read-only'):
            ref['nemitt_y'] = 4e-6
        with pytest.raises(AttributeError, match='read-only'):
            del ref['nemitt_y']
        with pytest.raises(AttributeError, match='read-only'):
            ref.beta_gamma_rel = 2.
        tw = xt.Table(dict(name=np.array(['c', 'end']), s=np.array([1., 1.6]),
                           betx=np.array([10., 10.]), bety=np.array([12., 12.]),
                           alfx=np.zeros(2), alfy=np.zeros(2),
                           x=np.zeros(2), y=np.zeros(2)))
        coll.assign_optics(twiss=tw, name='c', beta_gamma_rel=100.)
        ref = coll.reference_optics
        assert coll.optics_ready()
        assert ref.gemitt_x == pytest.approx(2.5e-8)
        assert ref.gemitt_y == pytest.approx(3.5e-8)
        assert len(ref) == 3
        jaw = coll.jaw_L if cls is xc.BlackAbsorber else coll.jaw_U
        coll.set_reference_emittance(nemitt_x=1e-5)
        assert (coll.jaw_L if cls is xc.BlackAbsorber else coll.jaw_U) == pytest.approx(2*jaw)
        assert coll.gap == pytest.approx(6.)
        copied = coll.copy()
        copied.set_reference_emittance(nemitt_x=2e-5)
        assert coll.reference_optics.nemitt_x == 1e-5
        assert copied.reference_optics.nemitt_x == 2e-5
        coll.set_reference_emittance(nemitt_y=None)
        assert not coll.optics_ready()
        assert coll.reference_optics.nemitt_y is None
        # Serialise the legacy emittance field; assigned optics rows are
        # deliberately not serialised. Check a pending partial round-trip.
        assert coll.to_dict()['emittance'] == [1e-5, None]
        pending = cls(length=0.6, gap=6.)
        pending.set_reference_emittance(nemitt_x=1e-5)
        restored = cls.from_dict(pending.to_dict())
        assert restored.reference_optics.nemitt_x == 1e-5
        assert restored.reference_optics.nemitt_y is None
    with pytest.warns(FutureWarning, match='nemitt_x'):
        assert coll.nemitt_x == 1e-5
    with pytest.warns(FutureWarning, match='nemitt_y'):
        coll.nemitt_y = 3.5e-6
    with pytest.warns(FutureWarning, match='gemitt_x'):
        assert coll.gemitt_x == pytest.approx(1e-7)
    with pytest.warns(FutureWarning, match='emittance'):
        coll.emittance = 2.5e-6
    assert coll.reference_optics.nemitt_y == 2.5e-6


@pytest.mark.xcother
def test_line_reference_emittance(make_line, monkeypatch):
    line = make_line()
    api = line.xcoll.collimators
    with warnings.catch_warnings():
        warnings.filterwarnings('error', category=FutureWarning)
        api.set_reference_emittance(nemitt_x=2.5e-6, nemitt_y=3.5e-6)
        assert len(api.reference_optics) == 0
        assert api.reference_optics.nemitt_x == 2.5e-6
        api.set_reference_emittance(names='coll', nemitt_x=4e-6)
        assert line['coll'].reference_optics.nemitt_x == 4e-6
        assert line['crystal'].reference_optics.nemitt_x == 2.5e-6
        assert len(api.reference_optics) == 0
        table = line.get_table()
        n = len(table.name)
        twiss = xt.Table(dict(name=table.name.copy(), s=table.s.copy(),
                              betx=np.full(n, 10.), bety=np.full(n, 12.),
                              alfx=np.zeros(n), alfy=np.zeros(n),
                              x=np.zeros(n), y=np.zeros(n)))
        def forbid_twiss(*args, **kwargs):
            pytest.fail('Precomputed optics or emittance updates must not compute Twiss')
        monkeypatch.setattr(xt.Line, 'twiss', forbid_twiss)
        api.assign_optics(twiss=twiss)
        ref = api.reference_optics
        assert 'nemitt_x' in ref._col_names
        assert np.all(ref.rows[ref.name == 'coll'].nemitt_x == 4e-6)
        assert ref.nemitt_y == 3.5e-6
        with pytest.raises(ValueError, match='read-only'):
            ref.nemitt_x[:] = 1e-6
        with pytest.raises(AttributeError, match='read-only'):
            ref['nemitt_y'] = 1e-6
        api.assign_optics(names='coll', twiss=twiss, nemitt_x=5e-6, nemitt_y=6e-6)
        assert line['coll'].reference_optics.nemitt_x == 5e-6
        assert line['coll'].reference_optics.nemitt_y == 6e-6
        assert line['crystal'].reference_optics.nemitt_y == 3.5e-6
        api.set_reference_emittance(nemitt_x=7e-6)
        assert api.reference_optics.nemitt_x == 7e-6
