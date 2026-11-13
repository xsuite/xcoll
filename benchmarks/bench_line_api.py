"""Guard the cost of copying an installed collimator (no tracking or optics)."""

import pytest
import xtrack as xt
import xcoll as xc


@pytest.mark.quick
@pytest.mark.benchmark(group='Line tools | element copy')
@pytest.mark.parametrize('installed', [False, True], ids=['standalone', 'installed'])
def bench_collimator_copy(benchmark, installed):
    coll = xc.BlackAbsorber(length=0.6, jaw=0.002)
    if installed:
        line = xt.Line(elements={
            'before': xt.Drift(length=2.),
            'coll': xt.Device(length=0.6),
            'after': xt.Drift(length=2.),
        })
        line.xcoll.collimators.install('coll', coll, need_apertures=False)
        coll = line.get('coll')

    copied = benchmark(coll.copy)
    assert copied is not coll
    assert copied.name == coll.name
    assert copied.jaw == coll.jaw
    if installed:
        assert coll._line is line
        assert copied._line is line
    else:
        assert '_line' not in copied.__dict__


@pytest.mark.quick
@pytest.mark.benchmark(group='Line tools | optics')
def bench_collimator_optics_at_positions(benchmark):
    import numpy as np

    # No tracking, compilation or Twiss calculation in the timed region.
    names = [f'coll{i}' for i in range(60)]
    line = xt.Line(elements={name: xc.BlackAbsorber(length=0.6) for name in names})
    table = line.get_table()
    n = len(table.name)
    twiss = xt.Table(dict(name=table.name.copy(), s=table.s.copy(),
                          betx=np.full(n, 10.), bety=np.full(n, 12.),
                          alfx=np.full(n, -0.5), alfy=np.zeros(n),
                          x=np.zeros(n), y=np.zeros(n)))
    optics = benchmark(line.xcoll.collimators.get_optics_at, names,
                       positions=['upstream', 'centre', 'downstream'], twiss=twiss)
    assert isinstance(optics, xt.Table)
    assert set(optics.name) == set(names)
    assert len(optics) == 3 * len(names)
    np.testing.assert_allclose(optics.s_rel.reshape(-1, 3),
                               np.tile([0., 0.3, 0.6], (len(names), 1)))
