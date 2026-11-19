import json
from pathlib import Path

import numpy as np
import pytest
import xtrack as xt
import xcoll as xc
from xcoll.initial_distribution import generate_delta_from_dispersion


@pytest.mark.parametrize('cls', [xc.BlackAbsorber, xc.BlackCrystal])
def test_collimator_geometric_optics(cls):
    line = xt.Line(elements={'coll': cls(length=0.1), 'drift': xt.Drift(length=1.)},
                   particle_ref=xt.Particles(p0c=2e9))
    coll = line.get('coll')
    bg = line.particle_ref.beta0[0] * line.particle_ref.gamma0[0]
    tt = line.get_table()
    tw = xt.Table(dict(name=tt.name, s=tt.s, betx=np.full(3, 10.), bety=np.full(3, 12.),
                       alfx=np.zeros(3), alfy=np.zeros(3), x=np.zeros(3), y=np.zeros(3)))
    norm = dict(nemitt_x=2e-6, nemitt_y=3e-6)
    geom = dict(gemitt_x=2e-6/bg, gemitt_y=3e-6/bg)
    line.xcoll.collimators.assign_optics(twiss=tw, **norm)
    coll.gap = 5
    jaw = coll.jaw
    line.xcoll.collimators.assign_optics(twiss=tw, **geom)
    np.testing.assert_allclose(coll.jaw, jaw)
    assert coll.reference_optics.nemitt_x == pytest.approx(2e-6)
    coll.assign_optics(reference_optics=coll.reference_optics,
                       gemitt_x=4e-6/bg, nemitt_y=3e-6)
    assert coll.reference_optics.nemitt_x == pytest.approx(4e-6)
    with pytest.raises(ValueError, match='Cannot provide both'):
        coll.assign_optics(reference_optics=coll.reference_optics, **norm, **geom)
    line.xcoll.collimators.set_reference_emittance(**geom)
    np.testing.assert_allclose(coll.jaw, jaw)
    coll.set_reference_emittance(gemitt_x=4e-6/bg)
    assert coll.reference_optics.nemitt_x == pytest.approx(4e-6)
    coll.set_reference_emittance(gemitt_x=None)
    assert coll.reference_optics.nemitt_x is None


def test_database_geometric_emittances():
    data = json.loads((Path(__file__).parent / 'data/colldb_lhc_run3_b1.json').read_text())
    ex, ey = data['emittance']['x'], data['emittance']['y']
    bg = 2.3
    original = xc.CollimatorDatabase.from_dict(data)
    geometric = xc.CollimatorDatabase.from_dict(data, gemitt_x=ex/bg,
                                                gemitt_y=ey/bg, beta_gamma_rel=bg)
    assert geometric.nemitt_x == pytest.approx(original.nemitt_x)
    assert geometric.nemitt_y == pytest.approx(original.nemitt_y)
    without_emittance = {k: v for k, v in data.items() if k != 'emittance'}
    mixed = xc.CollimatorDatabase.from_dict(without_emittance, gemitt_x=ex/bg,
                                            nemitt_y=ey, beta_gamma_rel=bg)
    assert mixed.nemitt_x == pytest.approx(ex)
    with pytest.raises(ValueError, match='beta_gamma_rel'):
        xc.CollimatorDatabase.from_dict(data, gemitt_x=ex/bg, gemitt_y=ey/bg)
    with pytest.raises(ValueError, match='Cannot provide both'):
        xc.CollimatorDatabase.from_dict(data, nemitt_x=ex, gemitt_x=ex/bg,
                                        nemitt_y=ey, beta_gamma_rel=bg)


def test_blowup_geometric_emittance():
    blowup = xc.BlowUp(plane='H')
    blowup.name = 'blowup'
    blowup.line = None
    tw = xt.Table(dict(name=np.array(['blowup']), gamx=np.array([0.1])))
    blowup.calibrate_by_emittance(2e-6, beta_gamma_rel=2.3, twiss=tw)
    expected = blowup.calibration
    blowup.calibrate_by_emittance(gemitt=2e-6/2.3, twiss=tw)
    assert blowup.calibration == pytest.approx(expected)
    with pytest.raises(ValueError, match='Cannot provide both'):
        blowup.calibrate_by_emittance(2e-6, gemitt=1e-6, twiss=tw)


def test_delta_geometric_emittances():
    line = xt.Line(elements={'coll': xt.Drift(length=0.1), 'rest': xt.Drift(length=1.)},
                   particle_ref=xt.Particles(p0c=2e9))
    line.build_tracker()
    tw = line.twiss4d(betx=3, bety=5, dx=0.5)
    bg = line.particle_ref.beta0[0] * line.particle_ref.gamma0[0]
    kwargs = dict(line=line, at_element='coll', plane='x', position_mm=0.002,
                  betatron_cut=2., twiss=tw)
    dn = generate_delta_from_dispersion(**kwargs, nemitt_x=2e-6, nemitt_y=3e-6)
    dg = generate_delta_from_dispersion(**kwargs, gemitt_x=2e-6/bg, gemitt_y=3e-6/bg)
    np.testing.assert_allclose(dn, dg, rtol=1e-13)
