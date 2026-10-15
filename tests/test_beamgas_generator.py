# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #
import json
import warnings

import numpy as np
import pandas as pd
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


COORDINATES = ('x', 'px', 'y', 'py', 'zeta', 'delta', 'weight', 's',
               'at_element', 'particle_id')


def _assert_same_particles(part, other):
    alive, alive_other = part.state > 0, other.state > 0
    assert np.sum(alive) == np.sum(alive_other) > 0
    for nn in COORDINATES:
        np.testing.assert_array_equal(
            getattr(part, nn)[alive], getattr(other, nn)[alive_other],
            err_msg=f'{nn} differs')


def test_generate_at_one_or_several_markers(ring):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1])
    bg = line.xcoll.beamgas
    names = bg.marker_names

    assert isinstance(bg.generate(at=names[1]), xt.Particles)
    assert list(bg.generate(at=[names[3], names[1]])) == [names[1], names[3]]
    particles, logs = bg.generate(return_log=True)
    assert list(particles) == list(logs) == names
    part, log = bg.generate(at=names[1], return_log=True)
    assert isinstance(log, xt.Table)
    assert line[names[1]].scatter_log is not None


@pytest.mark.parametrize('process', ['brems', 'coulomb'])
def test_generated_particles(ring, process):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1], process=process,
               coulomb_theta=(1e-6, 0.3))
    bg = line.xcoll.beamgas
    ref = line.particle_ref

    tab = line.get_table()
    for nn in bg.marker_names:
        part, log = bg.generate(at=nn, return_log=True)
        alive = part.state > 0
        assert np.sum(alive) == 200
        assert np.all(part.s[alive] == tab['s', nn])
        assert np.all(part.at_element[alive]
                      == line.element_names.index(nn))
        assert part.q0 == ref.q0 and part.mass0 == ref.mass0
        assert np.all(part.pdg_id[alive] == ref.pdg_id[0])
        # The weights carry the interaction rate of the section: exactly for
        # bremsstrahlung, on average for the importance-sampled Coulomb angles
        if process == 'brems':
            assert np.sum(part.weight[alive]) == pytest.approx(
                line[nn].interaction_rate, rel=1e-12)
        # One interaction per particle, on the gas species of the profile
        np.testing.assert_array_equal(np.sort(log.particle_id),
                                      np.sort(part.particle_id[alive]))
        assert set(log.gas) <= {'N', 'H'}
        assert np.all(log.name == nn)
        if process == 'coulomb':
            assert np.all((log.theta >= 1e-6) & (log.theta <= 0.3))
        else:
            assert np.all(log.photon_energy >= 1e6)


def test_generation_does_not_depend_on_the_other_markers(ring):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1])
    bg = line.xcoll.beamgas
    names = bg.marker_names

    single_job = bg.generate()
    for nn in reversed(names):
        _assert_same_particles(bg.generate(at=nn), single_job[nn])
    for nn, part in bg.generate(at=names[::-2]).items():
        _assert_same_particles(part, single_job[nn])


def test_replicas(ring):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1])
    bg = line.xcoll.beamgas
    nn = bg.marker_names[2]

    one = bg.generate(at=nn)
    first = bg.generate(at=nn, replica=0, n_replicas=3)
    replicas = [first] + [bg.generate(at=nn, replica=rr, n_replicas=3)
                          for rr in (1, 2)]

    # The same stream, weights divided by the number of replicas
    alive = one.state > 0
    np.testing.assert_array_equal(first.x[first.state > 0], one.x[alive])
    np.testing.assert_array_equal(
        first.weight[first.state > 0], one.weight[alive] * (1 / 3))
    # Independent samples, which together carry the rate once
    assert not np.array_equal(replicas[1].x[:10], replicas[2].x[:10])
    total = sum(np.sum(pp.weight[pp.state > 0]) for pp in replicas)
    assert total == pytest.approx(line[nn].interaction_rate, rel=1e-12)


@pytest.mark.parametrize('process', ['brems', 'coulomb'])
def test_generation_from_a_reloaded_line(ring, tmp_path, process):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1], process=process)
    names = line.xcoll.beamgas.marker_names
    expected = line.xcoll.beamgas.generate(at=names[1:4])

    line.to_json(tmp_path / 'line.json')
    line.env['ring'] = line
    line.env.to_json(tmp_path / 'env.json')
    for loaded in (xt.load(tmp_path / 'line.json'),
                   xt.load(tmp_path / 'env.json').lines['ring']):
        reloaded = loaded.xcoll.beamgas.generate(at=names[1:4])
        for nn, part in expected.items():
            _assert_same_particles(reloaded[nn], part)


def test_generate_checks_its_arguments(ring):
    line = _ring_with_markers()
    bg = line.xcoll.beamgas
    nn = bg.marker_names[1]
    with pytest.raises(ValueError, match='not configured'):
        bg.generate(at=nn)

    _configure(line, twiss=ring[1])
    with pytest.raises(ValueError, match='replica'):
        bg.generate(at=nn, replica=2, n_replicas=2)
    with pytest.raises(ValueError, match='n_replicas'):
        bg.generate(at=nn, n_replicas=0)

    line[nn].atomic_densities = {'N': 0.0}
    with pytest.raises(ValueError, match='gas densities are zero'):
        bg.generate(at=nn)

    line.set_particle_ref('electron', p0c=2e9)
    with pytest.raises(ValueError, match='configure the markers again'):
        bg.generate(at=nn)


def _fake_track(part, seed):
    """
    Stand-in for tracking: half of the particles are lost, on an aperture
    (state 0) or in a collimator (LOST_ON_MATERIAL = -330), and the first
    lost one produces a lost secondary, which carries its weight.
    """
    rng = np.random.default_rng(seed)
    alive = np.flatnonzero(part.state > 0)
    lost = alive[rng.random(len(alive)) < 0.5]
    part.state[lost] = rng.choice([0, -330], size=len(lost))
    part.at_turn[alive] = 10

    free = np.flatnonzero(part.particle_id < 0)[0]
    with part._bypass_linked_vars():
        for _, nn in xt.Particles.per_particle_vars:
            getattr(part, nn)[free] = getattr(part, nn)[lost[0]]
    part.particle_id[free] = np.max(part.particle_id) + 1
    part.parent_particle_id[free] = part.particle_id[lost[0]]


def _lost_contributions(part):
    primary = (part.parent_particle_id == part.particle_id) & (
        part.particle_id >= 0)
    lost = (part.state <= 0) & primary
    return part.weight[primary] * lost[primary]


def test_analyse_without_particles(ring):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1], process='coulomb',
               coulomb_theta=(1e-6, 0.3))
    bg = line.xcoll.beamgas
    result = bg.analyse()

    names = bg.marker_names
    rates = [line[nn].interaction_rate for nn in names]
    assert result.element_names == names
    assert result.rate_scattering == pytest.approx(np.sum(rates), rel=1e-12)
    assert result.lifetime_scattering == pytest.approx(
        BEAM['bunch_intensity'] / np.sum(rates), rel=1e-12)
    assert result.rate_above_theta_max == pytest.approx(
        sum(line[nn].rate_above_theta_max for nn in names), rel=1e-12)
    assert result.tracked is False
    assert result.rate_tracking is None and result.rate_generated is None
    assert set(result.local_rates._col_names) == {
        'name', 's', 'ds', 'mean_free_path', 'local_rate',
        'interaction_rate'}

    _configure(line, twiss=ring[1])
    assert bg.analyse().rate_above_theta_max is None


def test_analyse_merged_particles_and_samples_by_marker(ring):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1])
    bg = line.xcoll.beamgas
    samples = bg.generate()

    untracked = bg.analyse(xt.Particles.merge(list(samples.values())))
    assert untracked.tracked is False
    assert untracked.rate_tracking is None
    assert untracked.rate_generated == pytest.approx(
        untracked.rate_scattering, rel=1e-12)

    for ii, part in enumerate(samples.values()):
        _fake_track(part, seed=ii)
    expected = sum(np.sum(_lost_contributions(pp))
                   for pp in samples.values())

    by_marker = bg.analyse(samples)
    merged = bg.analyse(xt.Particles.merge(list(samples.values())[::-1]))
    singles = [bg.analyse(pp, at=nn) for nn, pp in samples.items()]

    for result in (by_marker, merged):
        assert result.tracked is True
        assert result.rate_tracking == pytest.approx(expected, rel=1e-12)
        assert result.lifetime_tracking == pytest.approx(
            BEAM['bunch_intensity'] / expected, rel=1e-12)
    assert sum(rr.rate_tracking for rr in singles) == pytest.approx(
        expected, rel=1e-12)

    # Per-marker counts and errors: from the samples of each marker only
    assert merged.rate_tracking_error is None
    assert 'sum_lost_weight' not in merged.local_rates._col_names
    variance = 0
    for single, nn in zip(singles, samples):
        contribution = _lost_contributions(samples[nn])
        n = len(contribution)
        error = np.sqrt(n * np.var(contribution, ddof=1))
        assert single.rate_tracking_error == pytest.approx(error, rel=1e-9)
        assert by_marker.local_rates['sum_lost_weight_error', nn] == (
            pytest.approx(error, rel=1e-9))
        for col in ('num_particles', 'sum_weight', 'num_lost_particles',
                    'sum_lost_weight'):
            assert single.local_rates[col, nn] == by_marker.local_rates[col, nn]
        variance += error**2
    assert by_marker.rate_tracking_error == pytest.approx(
        np.sqrt(variance), rel=1e-9)


def test_analyse_replicas(ring):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1])
    bg = line.xcoll.beamgas
    n_replicas = 3
    replicas = {nn: [] for nn in bg.marker_names}
    for rr in range(n_replicas):
        for nn, part in bg.generate(replica=rr,
                                    n_replicas=n_replicas).items():
            _fake_track(part, seed=rr)
            replicas[nn].append(part)

    result = bg.analyse(replicas)
    merged = bg.analyse(xt.Particles.merge(
        [pp for vv in replicas.values() for pp in vv]))
    assert result.rate_tracking == pytest.approx(merged.rate_tracking,
                                                 rel=1e-12)
    assert result.rate_generated == pytest.approx(result.rate_scattering,
                                                  rel=1e-12)
    np.testing.assert_array_equal(result.local_rates.n_replicas, n_replicas)
    # Independent samples: their variances add up
    variance = sum(len(cc) * np.var(cc, ddof=1)
                   for vv in replicas.values()
                   for cc in map(_lost_contributions, vv))
    assert result.rate_tracking_error == pytest.approx(np.sqrt(variance),
                                                       rel=1e-9)


@pytest.mark.parametrize('process', ['brems', 'coulomb'])
def test_analyse_with_the_records_of_the_interactions(ring, process):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1], process=process,
               coulomb_theta=(1e-6, 0.3))
    bg = line.xcoll.beamgas
    samples = {nn: [] for nn in bg.marker_names}
    logs = {nn: [] for nn in bg.marker_names}
    for rr in range(2):
        particles, records = bg.generate(replica=rr, n_replicas=2,
                                         return_log=True)
        for nn in particles:
            _fake_track(particles[nn], seed=rr)
            samples[nn].append(particles[nn])
            # Any table-like record works, e.g. one reloaded with pandas
            logs[nn].append(records[nn].to_pandas())

    with warnings.catch_warnings():
        # The fake tracking loses particles at all angles and energies: the
        # generation window looks truncated
        warnings.simplefilter('ignore')
        result = bg.analyse(samples, logs=logs)
    scan = result.cutoff_scan
    log = result.interaction_log
    assert len(log.particle_id) == 2 * 8 * 200
    assert set(log.replica) == {0, 1}
    # The first row of the scan is the whole sample
    assert scan.num_events.sum() == len(log.particle_id)
    assert scan.rate_tracking[0] == pytest.approx(result.rate_tracking,
                                                  rel=1e-12)
    assert scan.rate_tracking_error[0] == pytest.approx(
        result.rate_tracking_error, rel=1e-9)
    if process == 'coulomb':
        assert scan.cut[0] == pytest.approx(1e-6)
        assert result.rate_above_theta_max > 0
    else:
        assert scan.cut[0] == pytest.approx(1e6)
        assert result.rate_above_theta_max is None


def test_analyse_checks_its_arguments(ring):
    line = _ring_with_markers()
    bg = line.xcoll.beamgas
    nn = bg.marker_names[1]
    with pytest.raises(ValueError, match='not configured'):
        bg.analyse()

    _configure(line, twiss=ring[1])
    part, log = bg.generate(at=nn, return_log=True)
    other, other_log = bg.generate(at=nn, replica=1, n_replicas=2,
                                   return_log=True)
    with pytest.raises(ValueError, match='keys select the markers'):
        bg.analyse({nn: part}, at=nn)
    with pytest.raises(TypeError, match='xtrack.Particles'):
        bg.analyse({nn: [part, 'other']})
    with pytest.raises(ValueError, match='No particles'):
        bg.analyse({nn: []})
    with pytest.raises(ValueError, match='marker by marker'):
        bg.analyse(part, logs={nn: log})
    with pytest.raises(ValueError, match='same markers'):
        bg.analyse({nn: part}, logs={bg.marker_names[0]: log})
    with pytest.raises(ValueError, match='records of interactions'):
        bg.analyse({nn: [part, other]}, logs={nn: [log]})
    with pytest.raises(ValueError, match='does not match the particles'):
        bg.analyse({nn: part}, logs={nn: other_log})


def test_run_tracks_each_sample_from_its_marker(ring, monkeypatch):
    line = _ring_with_markers()
    _configure(line, twiss=ring[1])
    bg = line.xcoll.beamgas
    names = bg.marker_names

    calls = []

    def fake_track(particles, **kwargs):
        calls.append(kwargs)
        _fake_track(particles, seed=len(calls))

    monkeypatch.setattr(line, 'track', fake_track)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # see the test with the records
        result = bg.run(n_turns=7, at=names[1:4], n_replicas=2,
                        keep_particles=True)

    assert [cc['ele_start'] for cc in calls] == 2 * names[1:4]
    assert all(cc['ele_stop'] == cc['ele_start'] and cc['num_turns'] == 7
               for cc in calls)

    # The explicit steps give the same result
    samples = {nn: [] for nn in names[1:4]}
    for rr in range(2):
        for nn in names[1:4]:
            part = bg.generate(at=nn, replica=rr, n_replicas=2)
            _fake_track(part, seed=3 * rr + names[1:4].index(nn) + 1)
            samples[nn].append(part)
    explicit = bg.analyse(samples)
    assert result.element_names == explicit.element_names
    assert result.rate_tracking == pytest.approx(explicit.rate_tracking,
                                                 rel=1e-12)
    assert result.rate_tracking_error == pytest.approx(
        explicit.rate_tracking_error, rel=1e-12)
    assert result.cutoff_scan is not None
    assert len(result.interaction_log.particle_id) == 2 * 3 * 200

    assert list(result.particles_by_element) == names[1:4]
    assert np.all(result.lost_particles.state <= 0)


def test_run_tracks_the_particles(ring):
    line, tw = ring
    # No apertures in this ring, and tiny angles: nothing is lost
    _configure(line, twiss=tw, process='coulomb', coulomb_theta=(1e-6, 1e-5))
    result = line.xcoll.beamgas.run(n_turns=3, keep_particles=True)

    assert result.tracked is True
    assert result.rate_tracking == 0
    assert result.lifetime_tracking == np.inf
    assert result.lost_particles is None
    particles = result.particles
    assert np.all(particles.at_turn[particles.state > 0] == 3)


def test_run_checks_its_arguments(ring):
    with pytest.raises(ValueError, match='n_turns'):
        ring[0].xcoll.beamgas.run(n_turns=0)


def _assert_same_tracked_particles(part, other):
    """The same particles, in any order (matched by particle_id)."""
    def allocated_by_id(pp):
        index = np.flatnonzero(pp.particle_id >= 0)
        return index[np.argsort(pp.particle_id[index])]

    index, index_other = allocated_by_id(part), allocated_by_id(other)
    assert len(index) == len(index_other) > 0
    for nn in COORDINATES + ('state', 'at_turn', 'parent_particle_id'):
        np.testing.assert_array_equal(
            getattr(part, nn)[index], getattr(other, nn)[index_other],
            err_msg=f'{nn} differs')


def test_jobs_from_a_saved_line_merge_like_a_single_job(ring, tmp_path):
    """
    Independent jobs, each loading the saved line, generating one marker and
    replica, tracking and saving the particles and the record of their
    interactions, give the particles and the result of a single job.
    """
    line = _ring_with_markers()
    _configure(line, twiss=ring[1])
    line.to_json(tmp_path / 'line.json')
    bg = line.xcoll.beamgas
    names = bg.marker_names
    n_replicas = 2

    single_job = {nn: [] for nn in names}
    single_job_logs = {nn: [] for nn in names}
    for rr in range(n_replicas):
        particles, logs = bg.generate(replica=rr, n_replicas=n_replicas,
                                      return_log=True)
        for nn, part in particles.items():
            _fake_track(part, seed=rr)
            single_job[nn].append(part)
            single_job_logs[nn].append(logs[nn])

    # The jobs, in reverse order, exchanging only files
    jobs = [(nn, rr) for rr in range(n_replicas) for nn in names][::-1]
    for job_id, (nn, rr) in enumerate(jobs):
        job_line = xt.load(tmp_path / 'line.json')
        part, log = job_line.xcoll.beamgas.generate(
            at=nn, replica=rr, n_replicas=n_replicas, return_log=True)
        _fake_track(part, seed=rr)
        part.to_json(tmp_path / f'job.{job_id}.json')
        log.to_pandas().to_csv(tmp_path / f'job.{job_id}.csv', index=False)

    jobs_output = {nn: [None] * n_replicas for nn in names}
    jobs_logs = {nn: [None] * n_replicas for nn in names}
    for job_id, (nn, rr) in enumerate(jobs):
        with open(tmp_path / f'job.{job_id}.json') as fid:
            jobs_output[nn][rr] = xt.Particles.from_dict(json.load(fid))
        jobs_logs[nn][rr] = pd.read_csv(tmp_path / f'job.{job_id}.csv',
                                        float_precision='round_trip')

    for nn in names:
        for rr in range(n_replicas):
            _assert_same_tracked_particles(
                jobs_output[nn][rr], single_job[nn][rr])

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # see the test with the records
        expected = bg.analyse(single_job, logs=single_job_logs)
        result = bg.analyse(jobs_output, logs=jobs_logs)

    # Same particles, summed in another order (the loaded ones are
    # reorganized): equal to rounding
    assert result.rate_tracking == pytest.approx(
        expected.rate_tracking, rel=1e-12)
    assert result.rate_tracking_error == pytest.approx(
        expected.rate_tracking_error, rel=1e-9)
    for col in ('n_replicas', 'num_particles', 'num_lost_particles'):
        np.testing.assert_array_equal(
            result.local_rates[col], expected.local_rates[col])
    for col in ('sum_weight', 'sum_lost_weight'):
        np.testing.assert_allclose(
            result.local_rates[col], expected.local_rates[col], rtol=1e-12)
    for col in ('num_events', 'num_lost'):
        np.testing.assert_array_equal(
            result.cutoff_scan[col], expected.cutoff_scan[col])
    np.testing.assert_allclose(result.cutoff_scan.rate_tracking,
                               expected.cutoff_scan.rate_tracking, rtol=1e-12)
    for col in ('name', 'replica', 'particle_id', 'gas', 'photon_energy',
                'weight'):
        np.testing.assert_array_equal(
            result.interaction_log[col], expected.interaction_log[col])
    merged = bg.analyse(xt.Particles.merge(
        [pp for vv in jobs_output.values() for pp in vv]))
    assert merged.rate_tracking == pytest.approx(
        expected.rate_tracking, rel=1e-12)
