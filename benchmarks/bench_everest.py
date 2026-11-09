"""Standalone CPU performance benchmarks for Everest elements."""

import numpy as np
import pytest
import xobjects as xo
import xpart as xp
import xcoll as xc

N = 1_000_000
SEED_COORDS = 20261109
SEED_PARTICLES = 13579
ROUNDS = 5  # + 1 warm-up round => six tracking operations per case


def _record(benchmark, particles, element, label, start_px=None):
    """Diagnostics only; runs AFTER timing."""
    state = np.asarray(particles.state)
    lost = int(np.count_nonzero(state <= 0))
    benchmark.extra_info.update({
        "protocol": 2,
        "particles": N,
        "rounds": ROUNDS,
        "context": "ContextCpu-serial",
        "physics_case": label,
        "lost": lost,
        "lost_fraction": lost / N,
        "timed": "element.track(particles)",
        "coordinate_seed": SEED_COORDS,
        "particle_seed": SEED_PARTICLES,
    })
    if start_px is not None:
        ids = np.asarray(particles.particle_id, dtype=np.int64)
        valid = (state > 0) & (ids >= 0) & (ids < N)
        kick = np.asarray(particles.px)[valid] - start_px[ids[valid]]
        benchmark.extra_info["survivors_with_kick_gt_0.1urad"] = (
            float(np.mean(np.abs(kick) > 0.1e-6)) if kick.size else None
        )
        if label.startswith("crystal-channel"):
            benchmark.extra_info["channeling_proxy_fraction_of_survivors"] = (
                float(np.mean(kick > 0.5 * 149e-6)) if kick.size else None
            )
    return lost


@pytest.mark.quick
@pytest.mark.benchmark(group="Everest CPU | block")
def bench_block_cfc_1M(benchmark):
    """1 million protons at 450 GeV through 0.6 m CFC (bulk interactions)."""
    context = xo.ContextCpu()
    element = xc.EverestBlock(length=0.6, material=xc.materials.CarbonFibreComposite,
                              _context=context)
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=450e9, x=np.zeros(N), px=rng.normal(0, 3e-6, N),
        y=rng.normal(0, 1e-4, N), py=rng.normal(0, 3e-6, N),
        _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    element.track(particles0.copy())  # compile and warm up; excluded

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup, rounds=ROUNDS, warmup_rounds=1)
    assert _record(benchmark, result, element, "block-CFC", initial_px) > 0


@pytest.mark.full
@pytest.mark.benchmark(group="Everest CPU | block")
def bench_block_tungsten_1M(benchmark):
    """1 million 450 GeV protons through 0.6 m tungsten (thick-block limit)."""
    context = xo.ContextCpu()
    element = xc.EverestBlock(length=0.6, material=xc.materials.Tungsten,
                              _context=context)
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=450e9, x=np.zeros(N), px=rng.normal(0, 3e-6, N),
        y=rng.normal(0, 1e-4, N), py=rng.normal(0, 3e-6, N),
        _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    assert _record(benchmark, result, element, "block-W", initial_px) > 0


@pytest.mark.quick
@pytest.mark.benchmark(group="Everest CPU | collimator")
def bench_collimator_cfc_hit_1M(benchmark):
    """1 million 450 GeV protons inside the LEFT CFC jaw at entry (x > 1 mm)."""
    context = xo.ContextCpu()
    element = xc.EverestCollimator(
        length=0.6, jaw=[0.001, -0.001],
        material=xc.materials.CarbonFibreComposite, _context=context,
    )
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=450e9, x=rng.normal(1.08e-3, 5e-6, N),
        px=rng.normal(0, 3e-6, N), y=rng.normal(0, 1e-4, N),
        py=rng.normal(0, 3e-6, N), _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    assert np.all(np.asarray(particles0.x) > element.jaw_LU)
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    assert _record(benchmark, result, element, "collimator-CFC-hit", initial_px) > 0


@pytest.mark.full
@pytest.mark.benchmark(group="Everest CPU | collimator")
def bench_collimator_tungsten_hit_1M(benchmark):
    """1 million 450 GeV protons entering the left tungsten jaw at x=1.08 mm."""
    context = xo.ContextCpu()
    element = xc.EverestCollimator(
        length=0.6, jaw=[0.001, -0.001],
        material=xc.materials.Tungsten, _context=context,
    )
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=450e9, x=rng.normal(1.08e-3, 5e-6, N),
        px=rng.normal(0, 3e-6, N), y=rng.normal(0, 1e-4, N),
        py=rng.normal(0, 3e-6, N), _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    assert np.all(np.asarray(particles0.x) > element.jaw_LU)
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    assert _record(benchmark, result, element, "collimator-W-hit", initial_px) > 0


@pytest.mark.quick
@pytest.mark.benchmark(group="Everest CPU | collimator")
def bench_collimator_cfc_miss_1M(benchmark):
    """1 million 450 GeV protons inside the OPENING, far from either jaw."""
    context = xo.ContextCpu()
    element = xc.EverestCollimator(
        length=0.6, jaw=[0.001, -0.001],
        material=xc.materials.CarbonFibreComposite, _context=context,
    )
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=450e9, x=rng.normal(0, 1e-5, N),
        px=rng.normal(0, 3e-6, N), y=rng.normal(0, 1e-4, N),
        py=rng.normal(0, 3e-6, N), _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    assert np.max(np.abs(np.asarray(particles0.x))) < 0.001
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    assert _record(benchmark, result, element, "collimator-CFC-miss", initial_px) == 0


@pytest.mark.full
@pytest.mark.benchmark(group="Everest CPU | collimator")
def bench_collimator_cfc_grazing_1M(benchmark):
    """1 million protons near the left jaw; mixed shallow hits and near misses.

    More expensive than a deep hit is plausible: shallow particles can survive
    and scatter repeatedly, whereas particles entering deep may be absorbed.
    """
    context = xo.ContextCpu()
    element = xc.EverestCollimator(
        length=0.6, jaw=[0.001, -0.001],
        material=xc.materials.CarbonFibreComposite, _context=context,
    )
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=450e9, x=rng.normal(0.001, 2e-6, N),
        px=rng.normal(0, 6e-6, N), y=rng.normal(0, 1e-4, N),
        py=rng.normal(0, 3e-6, N), _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    initial_overjaw = float(np.mean(np.asarray(particles0.x) >= 0.001))
    assert 0.4 < initial_overjaw < 0.6
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    lost = _record(benchmark, result, element, "collimator-CFC-grazing", initial_px)
    assert lost > 0
    benchmark.extra_info["fraction_inside_jaw_at_entry"] = initial_overjaw


@pytest.mark.quick
@pytest.mark.benchmark(group="Everest CPU | crystal")
def bench_crystal_silicon_channeling_1M(benchmark):
    """1 million 400 GeV protons, 2 mm bent Si strip crystal (SPS-like reference)."""
    context = xo.ContextCpu()
    element = xc.EverestCrystal(
        length=0.002, material=xc.materials.Silicon,
        bending_angle=149e-6, width=0.002, height=0.05,
        side="+", jaw=0.001, miscut=0., lattice="strip",
        _context=context,
    )
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=400e9, x=rng.normal(1.5e-3, 75e-6, N),
        px=rng.normal(0, 5e-6, N), y=rng.normal(0, 1e-3, N),
        py=rng.normal(0, 5e-6, N), _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    _record(benchmark, result, element, "crystal-channel-Si", initial_px)
    proxy = benchmark.extra_info["channeling_proxy_fraction_of_survivors"]
    # A proxy is diagnostic; comparing it across branches is more meaningful
    # than imposing a hard threshold without a validated physics reference.
    assert proxy is not None


@pytest.mark.quick
@pytest.mark.benchmark(group="Everest CPU | crystal")
def bench_crystal_silicon_volume_reflection_1M(benchmark):
    """1 million 400 GeV protons, 2 mm bent Si strip crystal (SPS-like reference)."""
    context = xo.ContextCpu()
    element = xc.EverestCrystal(
        length=0.002, material=xc.materials.Silicon,
        bending_angle=149e-6, width=0.002, height=0.05,
        side="+", jaw=0.001, miscut=0., lattice="strip",
        _context=context,
    )
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=400e9, x=rng.normal(1.5e-3, 75e-6, N),
        px=rng.normal(75e-6, 15e-6, N), y=rng.normal(0, 1e-3, N),
        py=rng.normal(0, 5e-6, N), _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    _record(benchmark, result, element, "crystal-vr-Si", initial_px)


@pytest.mark.quick
@pytest.mark.benchmark(group="Everest CPU | crystal")
def bench_crystal_silicon_amorphous_1M(benchmark):
    """1 million 400 GeV protons, 2 mm bent Si strip crystal (SPS-like reference)."""
    context = xo.ContextCpu()
    element = xc.EverestCrystal(
        length=0.002, material=xc.materials.Silicon,
        bending_angle=149e-6, width=0.002, height=0.05,
        side="+", jaw=0.001, miscut=0., lattice="strip",
        _context=context,
    )
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=400e9, x=rng.normal(1.5e-3, 75e-6, N),
        px=np.concat([rng.normal(180e-6, 5e-6, N//2), rng.normal(40e-6, 5e-6, N-N//2)]),
        y=rng.normal(0, 1e-3, N),
        py=rng.normal(0, 5e-6, N), _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    _record(benchmark, result, element, "crystal-am-Si", initial_px)


@pytest.mark.full
@pytest.mark.benchmark(group="Everest CPU | crystal")
def bench_crystal_silicon_miss_1M(benchmark):
    """1 million protons well clear of a +1 mm bent crystal jaw."""
    context = xo.ContextCpu()
    element = xc.EverestCrystal(
        length=0.002, material=xc.materials.Silicon,
        bending_angle=149e-6, width=0.002, height=0.05,
        side="+", jaw=0.001, miscut=0., lattice="strip",
        _context=context,
    )
    rng = np.random.default_rng(SEED_COORDS)
    particles0 = xp.Particles(
        p0c=400e9, x=rng.normal(0, 0.1e-3, N),
        px=rng.normal(0, 30e-6, N), y=rng.normal(0, 1e-3, N),
        py=rng.normal(0, 5e-6, N), _context=context,
    )
    particles0._init_random_number_generator(
        seeds=SEED_PARTICLES + np.arange(N, dtype=np.uint32))
    initial_px = np.asarray(particles0.px).copy()
    element.track(particles0.copy())

    def setup():
        element._acc_ionisation_loss = 0.0
        element._acc_ionisation_loss_sec = 0.0
        return (particles0.copy(),), {}

    def track(particles):
        element.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup,
                                rounds=ROUNDS, warmup_rounds=1)
    assert _record(benchmark, result, element, "crystal-miss-Si", initial_px) == 0
