# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import pytest
import numpy as np
from pathlib import Path

import xtrack as xt
import xpart as xp
import xcoll as xc

from xobjects.test_helpers import for_all_test_contexts


path = Path(__file__).parent / 'data'
particle_ref = xt.Particles('proton', p0c=6.8e12)


def _bdsim_supports_tip_taper():
    # Tapered tips need a BDSIM including bdsim-collaboration/bdsim#228, and the interface compiled against it
    from g4interface import XtrackInterface
    return getattr(XtrackInterface, 'supportsTipTaper', lambda: False)()


@pytest.mark.geant4
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
def test_reload_bdsim(test_context):
    num_part = 1000
    _capacity = num_part*4

    if xc.geant4.engine.is_running():
        xc.geant4.engine.stop()

    coll = xc.Geant4Collimator(length=0.6, jaw=0.001, material='Ti', _context=test_context)
    xc.geant4.engine.particle_ref = particle_ref
    part_init = xp.build_particles(x=np.random.normal(coll.jaw_L + 1e-4, 1.e-4, num_part),
                                   particle_ref=particle_ref, _capacity=_capacity)
    part = []
    for _ in range(3):
        part.append(part_init.copy())
        cwd = Path.cwd()
        xc.geant4.engine.start(elements=coll, seed=1993)
        temp_cwd = xc.geant4.engine._cwd
        assert cwd != temp_cwd
        assert xc.geant4.engine._already_started
        assert (temp_cwd / 'rpyc.log').exists()
        assert (temp_cwd / 'root.out').exists()
        assert (temp_cwd / 'root.err').exists()
        assert (temp_cwd / 'geant4.out').exists()
        assert (temp_cwd / 'geant4.err').exists()
        assert (temp_cwd / 'engine.out').exists()
        assert (temp_cwd / 'engine.err').exists()
        assert not Path('rpyc.log').exists()
        assert not Path('root.out').exists()
        assert not Path('root.err').exists()
        assert not Path('geant4.out').exists()
        assert not Path('geant4.err').exists()
        assert not Path('engine.out').exists()
        assert not Path('engine.err').exists()

        t_start = time.time()
        coll.track(part[-1])
        print(f"Time per track: {(time.time()-t_start)*1e3:.2f}ms for "
            + f"{num_part} protons through {coll.length:.2f}m")
        assert (part[-1].state == xc.headers.particle_states.LOST_WITHOUT_SPEC).sum() == 0 # No particles should be lost without specification
        assert (part[-1].state == xc.headers.particle_states.LOST_ON_MATERIAL).sum() > 0   # Some particles should have died in the collimator
        assert (part[-1].state == 1).sum() > 0                                             # Some particles should have survived

        xc.geant4.engine.stop(clean=True)
        assert cwd == Path.cwd()
        assert not temp_cwd.exists()
        assert xc.geant4.engine._cwd is None
        assert xc.geant4.engine._g4link is None
        assert xc.geant4.engine._already_started
        assert not Path('rpyc.log').exists()
        assert not Path('root.out').exists()
        assert not Path('root.err').exists()
        assert not Path('geant4.out').exists()
        assert not Path('geant4.err').exists()
        assert not Path('engine.out').exists()
        assert not Path('engine.err').exists()

    # Check that the particles are the same
    for i in range(1, 3):
        assert np.allclose(part[0].x, part[i].x, atol=1e-12)
        assert np.allclose(part[0].px, part[i].px, atol=1e-12)
        assert np.allclose(part[0].y, part[i].y, atol=1e-12)
        assert np.allclose(part[0].py, part[i].py, atol=1e-12)
        assert np.allclose(part[0].zeta, part[i].zeta, atol=1e-12)
        assert np.allclose(part[0].delta, part[i].delta, atol=1e-12)
        assert np.all(part[0].state == part[i].state)
        assert np.all(part[0].particle_id == part[i].particle_id)
        assert np.all(part[0].parent_particle_id == part[i].parent_particle_id)


@pytest.mark.geant4
@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Geant4 only on CPU
)
def test_black_absorbers(test_context):
    n_part = 10_000
    _capacity = n_part*8
    angles = [0,45,90]
    angles = [0]
    jaws = np.array([0.03, -0.02])
    co = np.array([-0.01, 0.01])
    L = 0.873

    if xc.geant4.engine.is_running():
        xc.geant4.engine.stop()

    g4_collimators = []
    ba_collimators = []
    for angle in angles:
        shift = co[0]*np.cos(angle) + co[1]*np.sin(angle)
        g4coll = xc.Geant4Collimator(length=L, angle=angle, jaw=jaws+shift,
                                     _context=test_context, material='Beryllium')
        g4_collimators.append(g4coll)
        bacoll = xc.BlackAbsorber(length=L, angle=angle, jaw=jaws+shift,
                                  _context=test_context)
        ba_collimators.append(bacoll)

    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=g4_collimators, seed=1993, _all_black=True)

    x = np.random.uniform(-0.1, 0.1, n_part)
    y = np.random.uniform(-0.1, 0.1, n_part)
    px = np.random.uniform(-1e-3, 1e-3, n_part)
    py = np.random.uniform(-1e-3, 1e-3, n_part)
    part_init = xp.build_particles(x=x, y=y, px=px, py=py, _context=test_context,
                              particle_ref=xc.geant4.engine.particle_ref,
                              _capacity=_capacity)
    part = part_init.copy()
    part_ba = part_init.copy()

    for coll in g4_collimators:
        coll.track(part)

    for coll in ba_collimators:
        coll.track(part_ba)

    part.sort(interleave_lost_particles=True)
    part_ba.sort(interleave_lost_particles=True)

    assert np.all(part.state <= 1)
    assert np.all(part_ba.state <= 1)

    mask = part.state==1
    mask_ba = part_ba.state==1
    assert mask.sum() == mask_ba.sum()
    assert np.all(part.particle_id[mask] == part_ba.particle_id[mask_ba])

    # Stop the Geant4 connection
    xc.geant4.engine.stop(clean=True)


@pytest.mark.geant4
def test_geant4_tip_taper():
    if not _bdsim_supports_tip_taper():
        pytest.skip("BDSIM interface built without tip taper support (bdsim-collaboration/bdsim#228)")
    # SuperKEKB-type jaw, as in its engineering drawing and in SAD: a 10 mm long
    # flat face at the beam, a 12 degree taper over the first 37 mm of depth (up
    # to 358 mm long), then 88 mm more of jaw at that length. The element is
    # 10 mm long; Geant4 tracks through the full jaw, centred on it. The jaws
    # are deliberately not centred on the beam, to exercise the re-centring
    # that gives both jaws the same depth.
    if xc.geant4.engine.is_running():
        xc.geant4.engine.stop()

    length = 0.010
    jaw = [0.0080, -0.0076]
    jaw_depth = 0.125
    taper_depth = 0.037
    n_part = 400
    px = 1e-3
    kwargs = dict(length=length, jaw=jaw, jaw_depth=jaw_depth, taper_depth=taper_depth,
                  material=xc.materials.Copper,
                  tip_material=xc.materials.Tungsten, tip_thickness=0.005)

    flat = xc.Geant4CollimatorTip(**kwargs)
    tapered = xc.Geant4CollimatorTip(**kwargs, taper_angle_deg=12)
    assert flat.taper_extension == 0
    assert np.isclose(length + 2*tapered.taper_extension, 0.35814, atol=1e-5)

    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=[flat, tapered], seed=1993)

    def track(coll, depth):
        x0 = jaw[0] + depth
        part = xp.build_particles(x=np.full(n_part, x0), px=np.full(n_part, px),
                                  y=np.zeros(n_part), py=np.zeros(n_part),
                                  particle_ref=xc.geant4.engine.particle_ref,
                                  _capacity=n_part*4)
        coll.track(part)
        mask = (part.state == 1) & (part.particle_id < n_part)
        # Surviving primaries come out at the end of the element, not of the
        # longer wedge: x0 + px*length up to multiple scattering. Being 174 mm
        # off in either drift would shift this by px*174 mm = 174 um.
        if mask.sum() > 0:
            assert np.median(np.abs(part.x[mask] - (x0 + px*length))) < 30e-6
        return int(mask.sum())

    # 2 mm deep: 10 mm of tungsten through the flat jaw, ~29 mm through the taper
    # 30 mm deep: 10 mm of copper through the flat jaw, ~292 mm through the taper
    # 60 mm deep: 10 mm of copper through the flat jaw, the full 358 mm beyond the taper
    flat_shallow = track(flat, 0.002)
    flat_deep = track(flat, 0.030)
    flat_plateau = track(flat, 0.060)
    tapered_shallow = track(tapered, 0.002)
    tapered_deep = track(tapered, 0.030)
    tapered_plateau = track(tapered, 0.060)
    print(f"Surviving primaries out of {n_part} (2, 30, 60 mm deep): flat {flat_shallow}, "
          f"{flat_deep}, {flat_plateau}; tapered {tapered_shallow}, {tapered_deep}, "
          f"{tapered_plateau}")
    assert flat_deep > 0.8*n_part             # short jaw: deep hits mostly get through
    assert flat_plateau > 0.8*n_part
    assert tapered_deep < 0.5*flat_deep       # the taper stops them
    assert tapered_deep < tapered_shallow     # and stops deep hits more than shallow ones
    assert tapered_shallow < flat_shallow
    assert tapered_plateau <= tapered_deep    # beyond the taper the jaw stays at full length
    assert tapered_plateau < 0.25*n_part

    xc.geant4.engine.stop(clean=True)


@pytest.mark.geant4
def test_geant4_tip_jaw_width():
    if not _bdsim_supports_tip_taper():
        pytest.skip("BDSIM interface built without tip taper support (bdsim-collaboration/bdsim#228)")
    # The jaws are jaw_width wide across their plane of motion (BDSIM's default is
    # 12 mm). Hits further off-centre than half of that miss the jaw in Geant4.
    # A 20 mm wide jaw is used, so that hits at 8 mm only stop if it is passed on.
    if xc.geant4.engine.is_running():
        xc.geant4.engine.stop()

    jaw = [0.008, -0.008]
    n_part = 200
    coll = xc.Geant4CollimatorTip(length=0.010, jaw=jaw, jaw_depth=0.125, taper_depth=0.037,
                                  taper_angle_deg=12, jaw_width=0.020, material=xc.materials.Copper,
                                  tip_material=xc.materials.Tungsten, tip_thickness=0.005)

    xc.geant4.engine.particle_ref = particle_ref
    xc.geant4.engine.start(elements=coll, seed=1993)

    def survivors(y):
        # 60 mm deep in the jaw: 358 mm of copper, if inside it
        part = xp.build_particles(x=np.full(n_part, jaw[0] + 0.060), y=np.full(n_part, y),
                                  px=np.zeros(n_part), py=np.zeros(n_part),
                                  particle_ref=xc.geant4.engine.particle_ref,
                                  _capacity=n_part*4)
        coll.track(part)
        return int(((part.state == 1) & (part.particle_id < n_part)).sum())

    inside = survivors(0.008)
    outside = survivors(0.012)
    print(f"Surviving primaries out of {n_part}: {inside} at y=8 mm, {outside} at y=12 mm")
    assert inside < 0.25*n_part
    assert outside == n_part

    xc.geant4.engine.stop(clean=True)


@pytest.mark.geant4
def test_geant4_tip_taper_unsupported(monkeypatch):
    # An interface built without tip taper support must refuse a tapered tip, rather
    # than silently tracking through a flat jaw in a 2 m wide box
    from g4interface import XtrackInterface
    monkeypatch.setattr(XtrackInterface, 'supportsTipTaper', staticmethod(lambda: False),
                        raising=False)
    if xc.geant4.engine.is_running():
        xc.geant4.engine.stop()

    kwargs = dict(length=0.010, jaw=0.008, material=xc.materials.Copper,
                  tip_material=xc.materials.Tungsten, tip_thickness=0.005)
    xc.geant4.engine.particle_ref = particle_ref
    for extra in [dict(taper_angle_deg=12, taper_depth=0.037), dict(jaw_depth=0.125),
                  dict(jaw_width=0.02)]:
        with pytest.raises(RuntimeError, match="bdsim#228"):
            xc.geant4.engine.start(elements=xc.Geant4CollimatorTip(**kwargs, **extra), seed=1993)
        assert not xc.geant4.engine.is_running()

    # A tip without any of these is fine
    xc.geant4.engine.start(elements=xc.Geant4CollimatorTip(**kwargs), seed=1993)
    assert xc.geant4.engine.is_running()
    xc.geant4.engine.stop(clean=True)


def test_bdsim_source_gates(monkeypatch):
    # The markers in the interface source select lines by BDSIM version and by what
    # BDSIM provides (tip_taper), without needing BDSIM to be installed
    from xcoll.scattering_routines.geant4.environment import Geant4Interface
    src = '\n'.join([
        "always",
        "since_177  // BDSIM >= 1.7.7.develop",
        "call(a);   // BDSIM >= 1.7.7.develop BDSIM without tip_taper",
        "call(a,    // BDSIM with tip_taper",
        "b);        // BDSIM with tip_taper",
    ])

    def filtered(version, taper):
        monkeypatch.setattr(Geant4Interface, 'bdsim_has_feature',
                            lambda self, feature: taper)
        lines, _ = xc.geant4.interface._filter_version_gated_lines(src, version)
        return [line.split()[0] for line in lines]

    assert filtered('1.7.6', False) == ['always']
    assert filtered('1.7.7.develop', False) == ['always', 'since_177', 'call(a);']
    # 1.8.0.develop before and after bdsim-collaboration/bdsim#228: same version, different features
    assert filtered('1.8.0.develop', False) == ['always', 'since_177', 'call(a);']
    assert filtered('1.8.0.develop', True) == ['always', 'since_177', 'call(a,', 'b);']
