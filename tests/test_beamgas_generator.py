# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import numpy as np
import pytest

import xtrack as xt
import xcoll as xc


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
