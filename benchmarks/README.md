# Benchmarks

We use [pytest-benchmark](https://pytest-benchmark.readthedocs.io/) to measure the performance of our code while developing, and to keep track of how it evolves over time. We also store physical observables, such that we can check that a change in performance does not come at the expense of the physics.

## Quick start

Assuming you have a working Xcoll installation (typically your local checkout in editable mode), install pytest-benchmark in the same environment:

```bash
pip install pytest-benchmark
cd benchmarks
pytest -m quick                # A quick check while developing
```

Alternatively, from the repository root, use `pip install -e '.[benchmarks]'` to install Xcoll and the benchmark dependencies.

Some useful commands (always from the `benchmarks` folder):

```bash
pytest                         # All available benchmarks
pytest -m "not fluka and not geant4 and not gpu"  # All ordinary CPU benchmarks
pytest -m full                 # Only the expensive CPU cases (not the full suite!)
pytest -m fluka                # FLUKA benchmarks
pytest -m geant4               # Geant4/BDSIM benchmarks
pytest -m gpu                  # GPU benchmarks
pytest -k crystal              # Only cases with "crystal" in their name
pytest --collect-only -q       # List the benchmark cases
```

The results are saved automatically in `results/`. You can simply run the benchmarks again after changing your code and compare the two runs:

```bash
python -m perf.results list
python -m perf.results compare 2 1   # Second most recent vs most recent
python -m perf.report                # Make a local dashboard (results/dashboard.html)
```

Open the HTML file in your browser. There is also an [online dashboard](https://xsuite.github.io/xcoll/) with the results of the automated runs.

**Note:** do not use pytest-xdist (`-n` or `--tx`) for benchmarks. We run the benchmark cases one at a time, to avoid interference between them. When measuring a small change, try to use the same machine, Python environment, compiler, and Xsuite versions. Close other heavy applications, and avoid changing the benchmark parameters between the two measurements.

## How it works

The benchmarks are ordinary pytest functions using the `benchmark` fixture. Pytest-benchmark does the timing and writes its normal JSON output; `conftest.py` adds Xsuite-specific metadata. The tools to compare, store, and visualise results live in `perf/` and do not need to know the physics of the benchmarks. This is the part we can reuse in the other Xsuite packages.

### Adding a benchmark

Add a `bench_*` function in a `bench_*.py` file. The discovery rules, marks, and default pytest-benchmark options are defined in `pytest.ini`. Use a meaningful name, and give related cases the same group (e.g. `@pytest.mark.benchmark(group='Everest | tracking')`). Use `quick` for cases that are useful while developing, `full` for more expensive cases, and the appropriate engine or GPU marks for optional functionality.

As an example, using prepared `line` and `particles` fixtures:

```python
import numpy as np
import pytest

@pytest.mark.quick
@pytest.mark.benchmark(group='Tracking')
def bench_tracking(benchmark, line, particles):
    line.build_tracker()
    line.track(particles.copy(), num_turns=1)  # Warm-up, including compilation

    def setup():
        return (particles.copy(),), {}

    def track(p):
        line.track(p, num_turns=1)
        return p

    result = benchmark.pedantic(track, setup=setup, rounds=5)
    benchmark.extra_info.update(
        protocol=1, particles=len(particles.x), turns=1,
        lost_fraction=float(np.mean(result.state <= 0)),
    )
```

Here we measure tracking, not copying particles or compiling the kernel. If the purpose of a benchmark is to measure construction or compilation, then of course that belongs inside the timed section! The example assumes a fully populated particle array without secondaries; otherwise see the primary-ID diagnostic in `perf/particle_metrics.py`.

A few things to keep in mind:

- Reset mutable state between rounds (particles, counters, RNG state where appropriate). Record the seeds for stochastic tracking.
- Record workload and physics diagnostics in `benchmark.extra_info`. Increase `protocol` when the workload or timed operation changes, so old measurements are not compared as if they were identical.
- Optional engines can skip when unavailable. If a configured engine fails to start or track, the benchmark should fail. Engine start-up and interface compilation are outside the tracking measurement.

### Defining compatibility and physics checks

`metrics.toml` is the project-specific configuration. The `comparison.workload_fields` list defines which entries in `extra_info` describe the workload (number of particles, turns, geometry, seeds, protocol, etc.). An explicit `extra_info['workload']` dictionary can also be used. Two measurements are compatible if their workload dictionaries match; for timing we additionally require the same recorded machine identity.

The `[observables.*]` sections define which physical quantities to compare and their tolerances. For instance:

```toml
[observables.lost_fraction]
source = 'lost_fraction'
absolute_tolerance = 0.01
binomial_sigmas = 4
trials_field = 'particles'
```

For a lost fraction we allow a statistical fluctuation (here up to four binomial standard deviations), or the absolute tolerance, whichever is larger. Other quantities can use absolute and relative tolerances, or be normalised by another field (e.g. losses per particle). See the existing `metrics.toml` for examples. The comparison tools do not need to know what a collimator or particle is.

The physical diagnostics are a screening of selected scalar observables, not a substitute for our physics tests. Missing quantities are not considered agreement. The statuses distinguish an actual change (`CHANGED`) from an incompatible workload (`INCOMPATIBLE`), missing observables (`UNVERIFIED`), no available reference (`NOT_CHECKED`), a first scheduled measurement (`BASELINE`), and agreement within tolerance (`OK`).

### Local results and comparisons

`conftest.py` adds provenance to the native pytest-benchmark output: Xsuite and external-engine versions, machine identity, testbed, and compiler information. The result itself contains the timing samples and the Git revision. Use `BENCHMARK_MACHINE_NAME` if you want to give your machine an explicit stable name.

From `benchmarks/`:

```bash
python -m perf.results list
python -m perf.results tag 1 before-change
# Make a change and run the benchmark again
python -m perf.results tag 1 after-change
python -m perf.results compare before-change after-change
python -m perf.results untag before-change
```

The numbers refer to the listing (1 is the newest). Tags are local aliases in `results/.labels`, not Git tags; the JSON files are not renamed. Older files are archived automatically before pytest-benchmark reaches its four-digit file-number limit. To archive manually, use `python -m perf.results archive --force`.

The `compare` command runs pytest-benchmark's timing comparison and then the physics comparison. You can also use the standard pytest-benchmark CLI directly, or run just the physics check on two JSON files:

```bash
python -m perf.compare_physics old.json new.json
```

The dashboard (`python -m perf.report`) is a self-contained HTML file. Filter by branch, timeframe, or benchmark name; click a benchmark for its absolute runtimes, or use the comparison view to select two sets of measurements. The heatmap reference is selectable:

- **Previous scheduled (default):** the latest earlier compatible scheduled measurement, even outside the selected timeframe. Manual runs do not replace this reference. For local results without scheduled history, use the previous compatible measurement.
- **First, average, median:** use the compatible measurements in the selected timeframe, separately for each workload and machine. Average and median operate on the individual run medians.

The `vs first` column and absolute-runtime plot show longer-term evolution.

## Automated benchmarks and publication

The workflows are in `.github/workflows/`. There are three types of measurement:

| Workflow | What it runs | Where |
|---|---|---|
| `xcoll-performance.yml` | Full ordinary CPU suite daily (03:17 UTC); manual quick or full | Dedicated CERN VM for main; GitHub-hosted for other refs |
| `xcoll-performance-fluka.yml` | FLUKA daily (05:37 UTC) or manual | CERN VM |
| `xcoll-performance-geant4.yml` | Geant4/BDSIM daily (08:57 UTC) or manual | CERN VM |
| `xcoll-performance-pr.yml` | Full CPU main vs PR, sequentially on the same machine | GitHub-hosted |

GitHub schedules in UTC by default, and can start late. The CERN jobs share the same concurrency group, with a queue, to ensure that two benchmarks never run on the VM at the same time. The runner is restricted to reviewed main code. See [RUNNERS.md](RUNNERS.md) for the installation and security considerations. In particular, do not run arbitrary PR code on a persistent CERN runner.

Each CERN run creates a disposable Conda environment, installs released Xsuite and the current Xcoll checkout, and prebuilds the selected serial Xsuite kernels. FLUKA and Geant4 also compile their interfaces. The compiler-worker count (`0` for sequential compilation, `2` for two workers) is independent of the serial tracking context. Scheduled runs use `0`; manual runs can choose. The CPU workflow pins Conda GCC/G++ 16.2. We save the compiler and dependency versions and remove the environment afterwards. Released dependencies can change between days.

For FLUKA the machine needs a working installation (including system `tree` and `gfortran`), with the repository Actions variables `XCOLL_FLUKA_ROOT`, `XCOLL_FLUKA_SERVER`, `XCOLL_FLUKAIO_LIB`, `XCOLL_LINEBUILDER_ROOT`, and optionally `XCOLL_FLUKA_FEDB` pointing to the relevant paths. Geant4/BDSIM is installed in its Conda environment. See the respective workflows for the exact package specifications.

For PRs targeting main, we run the PR base revision first, then the PR merge commit on the **same** GitHub-hosted machine, in separate environments with shared dependency constraints and fresh kernels. This gives a direct performance comparison, unlike comparing to the CERN VM. The job summary and artefact contain timing changes and physics comparisons against both the paired main result and the scheduled CPU reference. New commits cancel older runs; draft PRs and `release/` PRs are skipped. No PR results enter permanent history. Benchmark failures fail the job, while physics comparison flags are informational.

### What is stored, and what is the baseline?

There are two levels of storage:

- GitHub Actions artefacts contain the detailed pytest-benchmark JSON and the dependency manifests (90 days). This also applies to PRs.
- The `benchmark-results` branch contains permanent, compact JSON records under `runs/`, and compressed rolling physics snapshots under `baselines/<physical-host>/`. These are written by the separate publication jobs, not by the benchmark runner. Individual timing samples and unnecessary bulky diagnostics are omitted from the permanent records.

Only the scheduled full CPU, FLUKA, and Geant4 runs replace their respective physics snapshots. The first successful run establishes a reference, and the next is compared to it before updating it. Manual runs read the latest scheduled snapshot (quick CPU uses the full CPU reference), but do not modify it. The PR comparison reads the CPU snapshot too. The snapshots contain configured observables and workload information, not particle arrays. This is a rolling reference, not an immutable scientific validation target: an intentional physics change will also become the next scheduled reference.

The website is regenerated from **all** permanent records. The three scheduled profiles on one host and UTC day share one column (latest successful run per profile); manual runs remain separate. Contributing revisions and timestamps are retained. History writes and deployments have separate concurrency queues. `xcoll-performance-publish.yml` is shared by the external engines; `xcoll-performance-dashboard.yml` deploys the existing history without running benchmarks. Pages must use **GitHub Actions** as its publishing source.

## Porting this to another Xsuite package

Most of the tooling can be reused as is. The important parts to adapt are the following:

1. Copy `benchmarks/perf/`, `conftest.py`, `pytest.ini`, `metrics.toml`, and `tests/pytest.ini`. Add a `[benchmarks]` extra with pytest and pytest-benchmark to the package's `pyproject.toml`.
2. Replace the `bench_*.py` files with benchmarks relevant to the package. Define the workload and the physical observables in `metrics.toml`, and adapt the version collection in `conftest.py` if needed. Use the same metadata convention for any new dependencies.
3. Adapt `prebuild_cpu.py` and the workflow setup to what needs to be compiled for the package. Some packages may not need all three tracking kernels, or may need other dependencies. `setup_external.py` and `particle_metrics.py` are only needed where relevant.
4. Adapt the workflow repository guards, branch/profile names, runner group and labels, testbed and physical-host IDs, environment setup, and benchmark selection. Reuse the history, physics, and report modules. The Xcoll online grouping assumes the `full`, `fluka`, and `geant4` profiles; adapt `ONLINE_PROFILES` in `report.py` and the baseline profile policy in `rolling_physics.py` to the new package. Do not blindly enable the Xcoll external-engine schedules.
5. Enable Pages with GitHub Actions, give the trusted publisher write access to the history branch, and set up the runner and its access restrictions. The first successful publication creates `benchmark-results`. Test the manual workflow before enabling scheduled execution.

Finally, keep the infrastructure tests. These check the storage, metadata, comparison, dashboard, and workflow logic without needing an external physics engine. Adapt the tests that refer to Xcoll-specific cases or workflow names. To run them (PyYAML is needed for the workflow tests; Node.js is optional for the JavaScript policy test):

```bash
pip install pytest pyyaml
cd benchmarks/tests
python -m pytest
```

The separate `pytest.ini` in this folder makes pytest collect `test_*`, not `bench_*`. This tests the infrastructure, not the performance of Xcoll.
