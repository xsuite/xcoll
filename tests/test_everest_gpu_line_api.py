# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np
import pytest
import xobjects as xo
import xtrack as xt
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts


@pytest.mark.gpu
@pytest.mark.xcother
@for_all_test_contexts(excluding=('ContextCpu', 'ContextPyopencl'))
def test_reference_optics_copy_to_gpu(test_context):
    # No full lattice or Twiss calculation is needed to exercise the optics
    # state, its Python-side copy, and the GPU-resident jaw fields.
    line = xt.Line(elements={
        'before': xt.Drift(length=1.),
        'coll': xt.Device(length=0.6),
        'crystal': xt.Device(length=0.004),
        'after': xt.Drift(length=1.),
    }, particle_ref=xt.Particles(p0c=7e12))
    line.xcoll.collimators.install(
        ['coll', 'crystal'],
        [xc.EverestCollimator(length=0.6, gap=6.,
                            material=xc.materials.MolybdenumGraphite),
         xc.EverestCrystal(length=0.004, gap=6.,
                         material=xc.materials.Silicon, bending_angle=50e-6)],
        apertures=xt.LimitEllipse(a=0.04, b=0.04), need_apertures=True)
    table = line.get_table()
    n = len(table.name)
    twiss = xt.Table(dict(name=table.name.copy(), s=table.s.copy(),
                          betx=np.full(n, 10.), bety=np.full(n, 12.),
                          alfx=np.zeros(n), alfy=np.zeros(n),
                          x=np.zeros(n), y=np.zeros(n)))
    line.xcoll.collimators.assign_optics(twiss=twiss, nemitt_x=2.5e-6,
                                       nemitt_y=3.5e-6)

    # Both a standalone element copy and an installed line copy must preserve
    # independent reference optics and the context-specific Everest state.
    copied_line = line.copy(_context=test_context)
    for name in ('coll', 'crystal'):
        source = line.get(name)
        source_ref = source.reference_optics
        original_jaw = source.jaw_U if name == 'crystal' else source.jaw_L
        original_gap = source.gap
        for copied in (source.copy(_context=test_context), copied_line.get(name)):
            assert copied._context is test_context
            assert copied.optics_ready()
            assert copied.reference_optics is not source_ref
            for column in ('position', 's', 's_rel', 'betx', 'bety', 'alfx', 'alfy', 'x', 'y'):
                np.testing.assert_array_equal(copied.reference_optics[column],
                                              source_ref[column])
            for key in ('nemitt_x', 'nemitt_y', 'beta_gamma_rel'):
                assert copied.reference_optics[key] == source_ref[key]
            assert copied.gap == pytest.approx(original_gap)
            assert copied.jaw == pytest.approx(source.jaw)
            # Quadrupling the horizontal emittance doubles a gap-controlled
            # horizontal jaw displacement (the mock closed orbit is zero).
            copied.set_reference_emittance(nemitt_x=1e-5)
            assert copied.gap == pytest.approx(original_gap)
            jaw = copied.jaw_U if name == 'crystal' else copied.jaw_L
            assert jaw == pytest.approx(2*original_jaw)
            assert source_ref.nemitt_x == 2.5e-6
            assert (source.jaw_U if name == 'crystal' else source.jaw_L) == original_jaw
            if name == 'coll':
                assert copied.jaw_R == pytest.approx(2*source.jaw_R)
            # Conversely, a jaw-controlled setting must not move when the
            # emittance changes; its effective gap must change instead.
            copied.gap = None  # Clear the manually set gap before setting a jaw.
            copied.jaw = 0.001
            fixed_jaw = copied.jaw
            gap = copied.gap
            copied.set_reference_emittance(nemitt_x=4e-5)
            assert copied.jaw == pytest.approx(fixed_jaw)
            assert copied.gap == pytest.approx(np.asarray(gap)/2, abs=1e-6)


@pytest.mark.gpu
@pytest.mark.xcother
@for_all_test_contexts(excluding=('ContextCpu', 'ContextPyopencl'))
def test_generate_pencil_on_copied_gpu_line(test_context):
    # A small stable ring gives real Twiss and particle construction without
    # compiling/tracking a full LHC line. Use 4D optics (no RF needed).
    cpu = xo.ContextCpu()
    line = xt.Line(elements={
        'map': xt.LineSegmentMap(length=1., betx=10., bety=12., qx=0.31, qy=0.32),
        'coll': xt.Device(length=0.6),
        'after': xt.Drift(length=0.4),
    }, particle_ref=xt.Particles(p0c=450e9))
    line.xcoll.collimators.install(
        'coll', xc.EverestCollimator(length=0.6, gap=6.,
                                    material=xc.materials.CarbonFibreComposite),
        apertures=xt.LimitEllipse(a=0.04, b=0.04), need_apertures=True)
    line.build_tracker(_context=cpu)
    twiss = line.twiss(method='4d')
    line.xcoll.collimators.assign_optics(twiss=twiss, nemitt_x=2.5e-6,
                                       nemitt_y=3.5e-6)
    copied = line.copy(_context=test_context)
    copied.build_tracker(_context=test_context)
    gpu_coll = copied.get('coll')
    assert gpu_coll._context is test_context
    assert gpu_coll.reference_optics is not line.get('coll').reference_optics
    assert gpu_coll.optics_ready()

    n = 512
    impact = 2e-6
    spread = 1e-6
    rng_state = np.random.get_state()
    try:
        np.random.seed(12345)
        for plane, angle in [('x', 0.), ('y', 90.)]:
            for converging in (True, False):
                for current in (line, copied):
                    coll = current.get('coll')
                    coll.angle = angle
                    coll.gap = 6.
                    sign = 1 if converging else -1
                    coll.tilt_L = coll.gap_L*coll.divergence + sign*1e-4
                    coll.tilt_R = coll.gap_R*coll.divergence - sign*1e-4
                # Default GPU context on the copied line, explicit GPU output
                # on a CPU line, and explicit CPU output from a GPU line.
                for current, kwargs, expected in (
                    (copied, {}, test_context),
                    (line, {'_context': test_context}, test_context),
                    (copied, {'_context': cpu}, cpu),
                ):
                    coll = current.get('coll')
                    part = current.xcoll.collimators.generate_pencil(
                        'coll', n, twiss=twiss, method='4d', side='+-',
                        impact_parameter=impact, pencil_spread=spread, **kwargs)
                    assert part._context is expected
                    ctx2np = expected.nparray_from_context_array
                    index = current.element_names.index('coll')
                    np.testing.assert_array_equal(ctx2np(part.at_element),
                                                  np.full(n, index))
                    assert part.start_tracking_at_element == index
                    assert np.all(ctx2np(part.state) > 0)
                    if not converging:
                        # The pencil was generated at the downstream jaw and
                        # backtracked to the entrance by the generator.
                        xt.Drift(length=coll.length, _context=expected).track(part)
                    coord = ctx2np(getattr(part, plane))
                    jaw_l = coll.jaw_LU if converging else coll.jaw_LD
                    jaw_r = coll.jaw_RU if converging else coll.jaw_RD
                    left = coord > 0
                    assert np.count_nonzero(left) == n//2
                    depth = np.where(left, coord-jaw_l, jaw_r-coord)
                    assert np.all(depth >= impact-1e-9)
                    assert np.all(depth <= impact+spread+1e-9)
                    assert np.all(np.isfinite(ctx2np(part.px)))
                    assert np.all(np.isfinite(ctx2np(part.py)))
    finally:
        np.random.set_state(rng_state)
