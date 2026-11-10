"""CPU benchmarks for the Xcoll InteractionRecord.

Run from benchmarks/ with::

    pytest bench_interaction_record.py

The timed tracking operations exclude input preparation, particle copying,
compilation, and record reset. Logging is reset before every round, and
record capacities are checked to avoid silently benchmarking an overflowed
record. The same particle distribution is used for each recording mode.
"""

import numpy as np
import pytest
import xobjects as xo
import xpart as xp
import xcoll as xc


@pytest.mark.quick
@pytest.mark.benchmark(group="InteractionRecord | transparent tracking")
@pytest.mark.parametrize("mode", ["off", "impacts", "impacts-exits", "compact"],
                         ids=["off", "impacts", "impacts-and-exits", "selected-columns"])
def bench_interaction_record_transparent(benchmark, mode):
    """Isolate logging overhead from material scattering using a transparent jaw.

    100k protons intersect the left jaw; tracking is otherwise deterministic.
    All modes track exactly the same particles.
    """
    n_particles = 100_000
    context = xo.ContextCpu()
    coll = xc.TransparentCollimator(
        name="TCP", length=0.6, jaw=[0.001, -0.001], _context=context,
    )
    rng = np.random.default_rng(20261110)
    particles0 = xp.Particles(
        p0c=450e9,
        x=rng.normal(1.025e-3, 2e-6, n_particles),
        px=rng.normal(0, 2e-6, n_particles),
        y=np.zeros(n_particles),
        py=np.zeros(n_particles),
        _context=context,
    )
    assert np.all(np.asarray(particles0.x) > coll.jaw_LU)

    record = None
    capacity = 3 * n_particles + 100
    if mode != "off":
        record = xc.InteractionRecord(
            elements=[coll], names=["TCP"], capacity=capacity,
            record_impacts=True,
            record_exits=mode in ("impacts-exits", "compact"),
            record_scatterings=False,
            columns=(["particle_id_before", "x_before"]
                     if mode == "compact" else None),
        )

    # Compile with the correct logging configuration before measuring.
    coll.track(particles0.copy())
    if record is not None:
        record._index.num_recorded = 0

    def setup():
        if record is not None:
            record._index.num_recorded = 0
        return (particles0.copy(),), {}

    def track(particles):
        coll.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup, rounds=5, warmup_rounds=1)
    lost = int(np.count_nonzero(np.asarray(result.state) <= 0))
    assert lost == 0

    recorded = int(record._index.num_recorded) if record is not None else 0
    if record is not None:
        assert 0 < recorded < capacity, (
            f"InteractionRecord contains {recorded} events (capacity {capacity}). "
            "Increase capacity to avoid truncation."
        )
    benchmark.extra_info.update({
        "protocol": 2,
        "particles": n_particles,
        "physics_case": "transparent-jaw-record",
        "record_mode": mode,
        "record_columns": "selected" if mode == "compact" else "all",
        "record_capacity": capacity if record is not None else 0,
        "recorded_interactions": recorded,
        "lost_fraction": 0.0,
        "timed": "coll.track(particles)",
        "coordinate_seed": 20261110,
    })
    if record is not None:
        record.stop(elements=[coll])


@pytest.mark.full
@pytest.mark.benchmark(group="InteractionRecord | Everest tracking")
@pytest.mark.parametrize("mode", ["off", "impacts-exits", "all-events", "selected-columns"],
                         ids=["off", "impacts-and-exits", "all-events", "selected-columns"])
def bench_interaction_record_everest(benchmark, mode):
    """Measure recording overhead with real CFC multiple scattering.

    15k protons enter a 0.6 m CFC jaw, with reproducible particle RNG states.
    The full event mode also records scattering events. All cases share the
    same physical input distribution; event recording is the only variable.
    """
    n_particles = 15_000
    context = xo.ContextCpu()
    coll = xc.EverestCollimator(
        name="TCP", length=0.6, jaw=[0.001, -0.001],
        material=xc.materials.CarbonFibreComposite, _context=context,
    )
    rng = np.random.default_rng(20261110)
    particles0 = xp.Particles(
        p0c=450e9,
        x=rng.normal(1.025e-3, 8e-6, n_particles),
        px=rng.normal(0, 3e-6, n_particles),
        y=rng.normal(0, 100e-6, n_particles),
        py=rng.normal(0, 3e-6, n_particles),
        _context=context,
    )
    particles0._init_random_number_generator(
        seeds=314159 + np.arange(n_particles, dtype=np.uint32)
    )

    record = None
    capacity = 600_000
    if mode != "off":
        record = xc.InteractionRecord(
            elements=[coll], names=["TCP"], capacity=capacity,
            record_impacts=True,
            record_exits=True,
            record_scatterings=mode in ("all-events", "selected-columns"),
            columns=(["particle_id_before", "x_before"]
                     if mode == "selected-columns" else None),
        )

    coll.track(particles0.copy())  # Compile and prewarm the chosen code path.
    if record is not None:
        record._index.num_recorded = 0

    def setup():
        coll._acc_ionisation_loss = 0.0
        coll._acc_ionisation_loss_sec = 0.0
        if record is not None:
            record._index.num_recorded = 0
        return (particles0.copy(),), {}

    def track(particles):
        coll.track(particles)
        return particles

    result = benchmark.pedantic(track, setup=setup, rounds=5, warmup_rounds=1)
    lost = int(np.count_nonzero(np.asarray(result.state) <= 0))
    assert 0 < lost < n_particles

    recorded = int(record._index.num_recorded) if record is not None else 0
    if record is not None:
        assert 0 < recorded < capacity, (
            f"InteractionRecord contains {recorded} events (capacity {capacity}). "
            "Increase capacity to avoid truncation."
        )
    benchmark.extra_info.update({
        "protocol": 2,
        "particles": n_particles,
        "physics_case": "everest-cfc-jaw-record",
        "record_mode": mode,
        "record_columns": "selected" if mode == "selected-columns" else "all",
        "record_capacity": capacity if record is not None else 0,
        "recorded_interactions": recorded,
        "lost_fraction": lost / n_particles,
        "timed": "coll.track(particles)",
        "coordinate_seed": 20261110,
        "particle_seed": 314159,
    })
    if record is not None:
        record.stop(elements=[coll])


@pytest.mark.full
@pytest.mark.benchmark(group="InteractionRecord | data extraction")
@pytest.mark.parametrize("format", ["full-jaw", "full-lattice", "selected-jaw"],
                         ids=["all-columns-jaw", "all-columns-lattice", "selected-columns-jaw"])
def bench_interaction_record_to_pandas(benchmark, format):
    """Measure conversion of a populated record, excluding event generation.

    The record is prepared once. Only record.to_pandas() is timed.
    """
    n_particles = 50_000
    context = xo.ContextCpu()
    coll = xc.TransparentCollimator(
        name="TCP", length=0.6, jaw=[0.001, -0.001], _context=context,
    )
    rng = np.random.default_rng(20261110)
    particles = xp.Particles(
        p0c=450e9,
        x=rng.normal(1.025e-3, 2e-6, n_particles),
        px=rng.normal(0, 2e-6, n_particles),
        y=np.zeros(n_particles),
        py=np.zeros(n_particles),
        _context=context,
    )

    selected = format == "selected-jaw"
    columns = ["particle_id_before", "x_before"] if selected else None
    frame = "lattice" if format == "full-lattice" else "jaw"
    capacity = 3 * n_particles + 100
    record = xc.InteractionRecord(
        elements=[coll], names=["TCP"], capacity=capacity,
        columns=columns, record_impacts=True, record_exits=True,
        record_scatterings=False,
    )
    coll.track(particles)
    recorded = int(record._index.num_recorded)
    assert 0 < recorded < capacity

    frame_data = benchmark.pedantic(
        record.to_pandas, kwargs={"frame": frame}, rounds=5, warmup_rounds=1
    )
    assert len(frame_data) == recorded
    assert "interaction_type" in frame_data.columns

    benchmark.extra_info.update({
        "protocol": 2,
        "particles": n_particles,
        "physics_case": "transparent-jaw-record-export",
        "record_mode": "impacts-exits",
        "record_columns": "selected" if selected else "all",
        "record_frame": frame,
        "record_capacity": capacity,
        "recorded_interactions": recorded,
        "output_columns": len(frame_data.columns),
        "timed": "record.to_pandas(frame=...)",
        "coordinate_seed": 20261110,
    })
    record.stop(elements=[coll])
