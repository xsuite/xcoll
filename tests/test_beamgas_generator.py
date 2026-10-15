# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import numpy as np
import pytest
from scipy import constants as sc

import xtrack as xt
import xcoll as xc
from xcoll.beamgas.cross_sections import (BremsstrahlungCalculator,
                                          CoulombScatteringCalculator)


def _ring():
    """Four-cell electron ring of the beam-gas example, without markers."""
    env = xt.Environment()
    line = env.new_line(components=[
        env.new('mqf.1', xt.Quadrupole, length=0.3, k1=0.1),
        env.new('d1.1', xt.Drift, length=1),
        env.new('mb1.1', xt.Bend, length=3, angle=np.pi/2),
        env.new('d2.1', xt.Drift, length=1),
        env.new('mqd.1', xt.Quadrupole, length=0.3, k1=-0.7),
        env.new('d3.1', xt.Drift, length=1),
        env.new('mb2.1', xt.Bend, length=3, angle=np.pi/2),
        env.new('d4.1', xt.Drift, length=1),
        env.new('mqf.2', xt.Quadrupole, length=0.3, k1=0.1),
        env.new('d1.2', xt.Drift, length=1),
        env.new('mb1.2', xt.Bend, length=3, angle=np.pi/2),
        env.new('d2.2', xt.Drift, length=1),
        env.new('mqd.2', xt.Quadrupole, length=0.3, k1=-0.7),
        env.new('d3.2', xt.Drift, length=1),
        env.new('mb2.2', xt.Bend, length=3, angle=np.pi/2),
        env.new('d4.2', xt.Drift, length=1),
    ])
    line.set_particle_ref('electron', p0c=1e9)
    line.configure_bend_model(core='full', edge=None)
    return line


BEAM = dict(nemitt_x=1e-5, nemitt_y=1e-7, sigma_z=4e-3, sigma_delta=1e-3,
            bunch_intensity=4e9)
DENSITY_N = 2 * 1e-7 * 1e2 / (sc.Boltzmann * 293.15)  # N2 at 1e-7 mbar


def _ring_with_markers():
    """The ring with eight equally spaced markers, ring closed."""
    line = _ring()
    length = line.get_length()
    line.xcoll.beamgas.install_markers(at=list(np.arange(1, 8) * length / 8))
    return line


def _gas(line):
    """A uniform N profile and a linearly increasing H one."""
    s = np.linspace(0, line.get_length(), 5)
    return xt.Table({
        'name': np.array([f'p{ii}' for ii in range(len(s))]),
        's': s,
        'N': np.full(len(s), DENSITY_N),
        'H': 1e13 * (1 + s),
    })


@pytest.fixture(scope='module')
def ring():
    """
    The ring with markers and its 4d twiss. Building a tracker compiles its
    kernels, for each new line: the other tests configure fresh rings with
    the same layout using this twiss, and track in this line.
    """
    line = _ring_with_markers()
    return line, line.twiss(method='4d')


def _configure(line, **kwargs):
    config = dict(gas_density=_gas(line), process='brems',
                  brems_energy_cut=1e6, n_scattering_events=200, seed=1997,
                  method='4d', verbose=False, **BEAM)
    config.update(kwargs)
    line.xcoll.beamgas.configure(**config)


def _place_markers(line, positions):
    """Install markers by hand, as in the examples: env.place + insert."""
    env = line.env
    places = []
    for ii, at in enumerate(positions):
        env.elements[f'bg.{ii}'] = xc.BeamGasMarker()
        places.append(env.place(f'bg.{ii}', at=at))
    line.insert(places)
    return [f'bg.{ii}' for ii in range(len(positions))]


def test_facade_returns_a_generator_bound_to_the_line():
    line = _ring()
    bg = line.xcoll.beamgas
    assert isinstance(bg, xc.BeamGasGenerator)
    assert bg.line is line


def test_generator_needs_a_line():
    with pytest.raises(TypeError, match='xtrack.Line'):
        xc.BeamGasGenerator(None)


def test_marker_names_are_in_line_order():
    line = _ring()
    assert line.xcoll.beamgas.marker_names == []

    length = line.get_length()
    names = _place_markers(line, [length, 3.0, 0.5 * length])
    assert line.xcoll.beamgas.marker_names == [names[1], names[2], names[0]]


def test_repeated_marker_is_refused():
    line = _ring()
    env = line.env
    env.elements['bg'] = xc.BeamGasMarker()
    line.insert([env.place('bg', at=1.0), env.place('bg', at=2.0)])
    with pytest.raises(ValueError, match='more than once'):
        line.xcoll.beamgas.marker_names


# The same positions given to env.place directly and to install_markers
POSITIONS = {
    'number':         lambda env: dict(at=3.0),
    'element':        lambda env: dict(at='mb1.1'),
    'start':          lambda env: dict(at='mqf.2@start'),
    'end':            lambda env: dict(at='mqd.1@end'),
    'from':           lambda env: dict(at=0.5, from_='mb2.1'),
    'from_anchor':    lambda env: dict(at=-0.1, from_='mb2.1',
                                       from_anchor='end'),
    'from_at_anchor': lambda env: dict(at=0.2, from_='mb1.2@start'),
    'expression':     lambda env: dict(at='2 * s_bg'),
    'reference':      lambda env: dict(at=3 * env.ref['s_bg']),
}


@pytest.mark.parametrize('position', POSITIONS.values(), ids=POSITIONS.keys())
def test_install_markers_places_like_env_place(position):
    manual = _ring()
    env = manual.env
    env['s_bg'] = 2.6
    env.elements['bg'] = xc.BeamGasMarker()
    manual.insert([env.place('bg', **position(env))])

    line = _ring()
    line.env['s_bg'] = 2.6
    names = line.xcoll.beamgas.install_markers(
        close_ring=False, **position(line.env))

    assert names == ['beamgas.0']
    assert isinstance(line['beamgas.0'], xc.BeamGasMarker)
    tab, tab_manual = line.get_table(), manual.get_table()
    assert tab['s', 'beamgas.0'] == tab_manual['s', 'bg']
    # The rest of the line is cut in the same way
    np.testing.assert_array_equal(tab.s, tab_manual.s)
    np.testing.assert_array_equal(tab.element_type, tab_manual.element_type)


def test_install_markers_closes_the_ring():
    line = _ring()
    bg = line.xcoll.beamgas

    names = bg.install_markers(at=[1.0, 5.0])
    assert names == ['beamgas.0', 'beamgas.1', 'beamgas.2']
    assert bg.marker_names == names
    assert line.get_table()['s', 'beamgas.2'] == line.get_length()

    # The ring is already closed: no other marker at the end
    assert bg.install_markers(at=8.0) == ['beamgas.3']
    assert bg.marker_names[-1] == 'beamgas.2'


def test_install_markers_with_a_marker_at_the_end():
    line = _ring()
    names = line.xcoll.beamgas.install_markers(at=[1.0, line.get_length()])
    assert names == ['beamgas.0', 'beamgas.1']


def test_install_markers_without_closing_the_ring():
    line = _ring()
    names = line.xcoll.beamgas.install_markers(
        at=[1.0, 5.0], close_ring=False)
    assert names == ['beamgas.0', 'beamgas.1']
    assert line.get_table()['s', 'beamgas.1'] < line.get_length()


def test_install_markers_names():
    line = _ring()
    bg = line.xcoll.beamgas
    line.env.elements['beamgas.1'] = xt.Marker()

    assert bg.install_markers(
        at=[1.0, 2.0], close_ring=False) == ['beamgas.0', 'beamgas.2']
    assert bg.install_markers(
        at=3.0, names='my_bg', close_ring=False) == ['my_bg']
    assert bg.install_markers(
        at=4.0, prefix='bg', close_ring=False) == ['bg.0']
    with pytest.raises(ValueError, match='already used'):
        bg.install_markers(at=5.0, names='my_bg')
    with pytest.raises(ValueError, match='repeated'):
        bg.install_markers(at=[5.0, 6.0], names=['a', 'a'])


def test_install_markers_checks_its_arguments():
    line = _ring()
    bg = line.xcoll.beamgas
    with pytest.raises(ValueError, match='`from_` has 1 entries'):
        bg.install_markers(at=[1.0, 2.0], from_=['mb1.1'])
    with pytest.raises(ValueError, match='`names` has 1 entries'):
        bg.install_markers(at=[1.0, 2.0], names=['a'])
    with pytest.raises(ValueError, match='at least one'):
        bg.install_markers(at=[])
    assert bg.marker_names == []


@pytest.mark.parametrize('process', ['brems', 'coulomb'])
def test_configure_stores_the_configuration_on_the_markers(ring, process):
    _, tw = ring
    line = _ring_with_markers()
    gas = _gas(line)
    _configure(line, twiss=tw, gas_density=gas, process=process,
               coulomb_theta=(1e-6, 0.3))
    bg = line.xcoll.beamgas
    names = bg.marker_names

    ref = line.particle_ref
    p0c = ref.p0c[0]
    if process == 'brems':
        xsec = {Z: BremsstrahlungCalculator(Z, p0c, energy_cut=1e6)
                .compute_xsec() for Z in (7, 1)}
    else:
        xsec = {Z: CoulombScatteringCalculator(Z, p0c, q0=-1,
                                               theta_lim=(1e-6, 0.3))
                .compute_xsec() for Z in (7, 1)}
    f_rev = ref.beta0[0] * sc.c / line.get_length()
    tab = line.get_table()

    s_before = 0.0
    for nn in names:
        marker = line[nn]
        s = tab['s', nn]
        assert marker.process == process
        assert marker.p0c == p0c
        assert marker.ds == pytest.approx(s - s_before)
        for name, col in [('betx', 'betx'), ('dx', 'dx'), ('x_co', 'x'),
                          ('delta_co', 'delta')]:
            assert getattr(marker, name) == pytest.approx(
                tw[col, nn], rel=1e-12, abs=1e-15)
        densities = marker.atomic_densities
        assert densities['N'] == pytest.approx(DENSITY_N)
        assert densities['H'] == pytest.approx(1e13 * (1 + s))
        assert marker.mean_free_path == pytest.approx(
            1 / (densities['N'] * xsec[7] + densities['H'] * xsec[1]))
        # The profiles are linear: the integrals over the section are exact
        integral_n = DENSITY_N * (s - s_before)
        integral_h = 1e13 * ((s - s_before) + (s**2 - s_before**2) / 2)
        assert marker.interaction_rate == pytest.approx(
            BEAM['bunch_intensity'] * f_rev
            * (integral_n * xsec[7] + integral_h * xsec[1]), rel=1e-10)
        assert marker.bunch_intensity == BEAM['bunch_intensity']
        assert marker.n_scattering_events == 200
        assert marker.seed == 1997
        if process == 'brems':
            assert marker.brems_energy_cut == 1e6
            assert marker.rate_above_theta_max == 0
        else:
            assert (marker.coulomb_theta_min, marker.coulomb_theta_max) == (
                1e-6, 0.3)
            assert marker.rate_above_theta_max > 0
        s_before = s


def test_configure_does_not_modify_the_gas_density(ring):
    line = _ring_with_markers()
    gas = _gas(line)
    h_before = gas.H.copy()
    _configure(line, twiss=ring[1], gas_density=gas)
    np.testing.assert_array_equal(gas.H, h_before)


def test_configure_stores_a_drawn_seed(ring):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1], seed=None)
    seeds = {line[nn].seed for nn in line.xcoll.beamgas.marker_names}
    assert len(seeds) == 1
    assert seeds.pop() > 0

    with pytest.raises(ValueError, match='non-negative'):
        _configure(line, twiss=ring[1], seed=-1)
    with pytest.raises(TypeError, match='integer'):
        _configure(line, twiss=ring[1], seed=1.5)


def test_configure_a_subset_of_the_markers(ring):
    line, tw = _ring_with_markers(), ring[1]
    _configure(line, twiss=tw)
    names = line.xcoll.beamgas.marker_names
    rates = {nn: line[nn].interaction_rate for nn in names}

    # One section integrated on its own gives the same rate
    for nn in names:
        line[nn].interaction_rate = 0
        _configure(line, twiss=tw, at=nn)
        assert line[nn].interaction_rate == pytest.approx(rates[nn],
                                                          rel=1e-12)

    # Only the selected markers are configured: the rate goes as N
    _configure(line, twiss=tw, at=names[2:4],
               bunch_intensity=2 * BEAM['bunch_intensity'])
    for nn in names:
        factor = 2 if nn in names[2:4] else 1
        assert line[nn].interaction_rate == pytest.approx(
            factor * rates[nn], rel=1e-12)


def test_configure_computes_the_twiss(ring):
    line, tw = ring
    _configure(line)
    computed = {nn: line[nn].to_dict()
                for nn in line.xcoll.beamgas.marker_names}
    _configure(line, twiss=tw)
    for nn, dct in computed.items():
        assert line[nn].to_dict().keys() == dct.keys()
        for kk, vv in dct.items():
            np.testing.assert_array_equal(line[nn].to_dict()[kk], vv)


def test_configure_requires_the_last_marker_at_the_end_of_the_line():
    line = _ring()
    line.xcoll.beamgas.install_markers(at=[1.0, 5.0], close_ring=False)
    with pytest.raises(ValueError, match='must sit at the end of the line'):
        _configure(line)


def test_configure_checks_its_arguments(ring):
    line = ring[0]
    with pytest.raises(ValueError, match='`sigma_z` is required'):
        _configure(line, sigma_z=None)
    with pytest.raises(ValueError, match='not both'):
        _configure(line, gemitt_x=1e-9, gemitt_y=1e-11)
    with pytest.raises(TypeError, match=r'xt\.Table'):
        _configure(line, gas_density={'s': [0.0]})
    with pytest.raises(ValueError, match='process'):
        _configure(line, process='other')
    with pytest.raises(ValueError, match="'mb1.1' is not present"):
        _configure(line, at='mb1.1')
    with pytest.raises(TypeError, match='not a BeamGasMarker'):
        _configure(line, at='mqf.1')
    for at in (42, [42]):
        with pytest.raises(TypeError, match='selected by name'):
            _configure(line, at=at)

    protons = _ring_with_markers()
    protons.set_particle_ref('proton', p0c=1e9)
    with pytest.raises(ValueError, match='electron and positron'):
        _configure(protons, twiss=ring[1])

    no_markers = _ring()
    with pytest.raises(ValueError, match='install_markers'):
        _configure(no_markers)
