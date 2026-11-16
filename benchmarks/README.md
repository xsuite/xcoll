# Benchmarks

We use [pytest-benchmark](https://pytest-benchmark.readthedocs.io/) to measure the performance of our code while developing, and to keep track of how it evolves over time. We also store physical observables, such that we can check that a change in performance does not come at the expense of the physics. In addition to these benchmark diagnostics, we maintain an accepted Everest physics reference for more detailed regression testing (see [Physics validation](#physics-validation)).

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

- Write generated files (plots, input files, logs, etc.) to `tmp_path`, `results/`, or another Git-ignored location. Do not create or modify non-ignored files during a benchmark: the Git provenance safeguard will reject the run.
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

Benchmark sessions that save results check the Git checkout at the beginning and immediately before saving. Uncommitted changes are allowed, but the commit, branch, staged and unstaged changes, and non-ignored untracked files (including their contents) must remain unchanged throughout the run. Do not edit, commit, or switch branches while benchmarks are running. Keep generated output in `results/`, another Git-ignored directory, or a temporary directory outside the checkout. An unchanged dirty tree is allowed for local development; commit hashes alone do not identify its changes. The check compares two snapshots, so changes made and reverted in between are not detected. Use an isolated checkout for reproducible runs.

The trusted history publisher also checks that the JSON's commit matches the commit recorded by the workflow before writing a permanent history record.

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

Heatmap colours show changes in median runtime on a fixed logarithmic ratio scale: teal for faster, rust/red for slower, and neutral within ±3%. The strongest colours correspond to a fourfold speed-up (−75%) or slowdown (+300%); a white cross marks cells at or beyond these limits. Hover for the exact percentage, including changes outside the colour scale. The red corner separately flags a physics change.

The `vs first` column and absolute-runtime plot show longer-term evolution.

## Automated benchmarks and publication

The workflows are in `.github/workflows/`. Scheduled benchmarks, PR benchmarks, and physics validation are separate workflows:

| Workflow | What it runs | Where |
|---|---|---|
| `xcoll-performance.yml` | Full ordinary CPU suite daily (03:17 UTC); manual quick or full | Dedicated CERN VM for main; GitHub-hosted for other refs |
| `xcoll-performance-fluka.yml` | FLUKA daily (05:37 UTC) or manual | CERN VM |
| `xcoll-performance-geant4.yml` | Geant4/BDSIM daily (08:57 UTC) or manual | CERN VM |
| `xcoll-performance-pr.yml` | Full CPU PR base vs PR merge, sequentially on the same machine | GitHub-hosted |
| `xcoll-physics-pr.yml` | Everest physics against the committed reference, automatically on PRs or manually by PR number | GitHub-hosted |
| `xcoll-performance-external-pr.yml` | Paired FLUKA or Geant4 PR benchmarks, manual only, administrator-authorised | CERN VM |

GitHub schedules in UTC by default, and can start late. The CERN jobs share the same concurrency group, with a queue, to ensure that two benchmarks never run on the VM at the same time. The scheduled workflows run reviewed main code. Manual external PR benchmarks use the same CERN VM, but only after the explicit administrator approval described below. See [RUNNERS.md](RUNNERS.md) for the installation and security considerations, and [PR_WORKFLOWS.md](PR_WORKFLOWS.md) for the additional PR workflow configuration. In particular, do not run arbitrary PR code on a persistent CERN runner.

Each CERN run creates a disposable Conda environment, installs released Xsuite and the current Xcoll checkout, and prebuilds the selected serial Xsuite kernels. FLUKA and Geant4 also compile their interfaces. The compiler-worker count (`0` for sequential compilation, `2` for two workers) is independent of the serial tracking context. Scheduled runs use `0`; manual runs can choose. The CPU workflow pins Conda GCC/G++ 16.2. We save the compiler and dependency versions and remove the environment afterwards. Released dependencies can change between days.

For FLUKA the machine needs a working installation (including system `tree` and `gfortran`), with the repository Actions variables `XCOLL_FLUKA_ROOT`, `XCOLL_FLUKA_SERVER`, `XCOLL_FLUKAIO_LIB`, `XCOLL_LINEBUILDER_ROOT`, and optionally `XCOLL_FLUKA_FEDB` pointing to the relevant paths. Geant4/BDSIM is installed in its Conda environment. See the respective workflows for the exact package specifications.

For the automatic CPU benchmarks on PRs targeting main, we run the PR base revision first, then the PR merge commit on the **same** GitHub-hosted machine, in separate environments with shared dependency constraints and fresh kernels. This gives a direct performance comparison, unlike comparing to the CERN VM. The job summary and artefact contain timing changes and physics comparisons against both the paired main result and the scheduled CPU reference. New commits cancel older runs; draft PRs and `release/` PRs are skipped. Benchmark failures fail the job, while physics comparison flags are informational.

Successful paired comparisons are also saved under `pr-runs/<PR number>/` on `benchmark-results` and shown in the **Pull Requests** tab of the online dashboard. Select a PR and a workflow run to inspect its saved comparison, including the workload parameters and physical observables. The results are frozen as measured: revisiting an old PR does not compare it again against today's physics baseline. The PR measurement jobs use read-only repository permissions; separate, trusted `workflow_run` publishers validate the artefacts, write the compact records, and update the website. The raw timing samples remain available only in the 90-day Actions artefact.

### Physics and external-engine PR workflows

The Everest reference-validation workflow is independent of the CPU performance workflow. It runs automatically for eligible PRs targeting main, and can also be dispatched manually from main with a PR number:

```bash
gh workflow run xcoll-physics-pr.yml --repo xsuite/xcoll --ref main -f pr_number=123
```

A manual dispatch tests the existing PR merge commit, without requiring the contributor's branch or GitHub's merge ref to be updated to the latest main. Automatic runs retain the merge revision from the PR event. Both use a separately pinned checkout of main for the accepted reference and validation tools, and compile kernels from the candidate checkout. No reference generation is performed in PR CI. A manual physics dispatch does not restart the CPU performance workflow; the physics jobs have their own concurrency group.

The FLUKA and Geant4 PR benchmarks are manual only. Once the administrator-only workflow has been configured and enabled (see [PR_WORKFLOWS.md](PR_WORKFLOWS.md)), a repository administrator can dispatch a run from main with the PR number, the engine, and the full reviewed head SHA:

```bash
gh workflow run xcoll-performance-external-pr.yml --repo xsuite/xcoll --ref main \
    -f pr_number=123 -f engine=fluka \
    -f approved_head_sha=REPLACE_WITH_REVIEWED_40_CHARACTER_SHA
```

Use `engine=geant4` for Geant4/BDSIM. These workflows require a merge against current main, the exact reviewed PR head, and the configured approval gate before scheduling on CERN. Both revisions are benchmarked sequentially in separate environments on the existing VM, with the same dependency configuration and freshly compiled kernels and interfaces. The CERN jobs share the scheduled benchmark queue; no additional runner is needed. Administrator approval is not a sandbox: approve only code that you trust to execute on the persistent VM.

The supplementary publisher automatically saves physics validation (including failed or incomplete runs) and successful external-engine PR comparisons in the same permanent PR history, and redeploys the dashboard. No manual history update is needed. The **Pull Requests** view distinguishes CPU, FLUKA, Geant4, and Everest physics runs. Results for identical base, head, and merge SHAs are linked; results from different tested revisions are labelled accordingly. The physics panel shows the per-case outcome, duration, and available screening values and limits, together with the reference-source and validation-tooling SHAs. See [PR_WORKFLOWS.md](PR_WORKFLOWS.md) for details, including the separate status posted for a manual physics run when its tested merge is still current.

### What is stored, and what is the baseline?

For performance benchmarks and PR results there are two levels of storage:

- GitHub Actions artefacts contain the detailed pytest-benchmark JSON and the dependency manifests (90 days). This also applies to PR benchmarks. The separate Everest physics workflow uploads JUnit reports, including the captured screening results.
- The `benchmark-results` branch contains permanent, compact JSON records under `runs/`, frozen PR comparisons and Everest validation records under `pr-runs/<PR number>/`, and compressed rolling physics snapshots under `baselines/<physical-host>/`. These are written by the separate publication jobs, not by the benchmark runner. Individual timing samples and unnecessary bulky diagnostics are omitted from the permanent records.

Only the scheduled full CPU, FLUKA, and Geant4 runs replace their respective physics snapshots. The first successful run establishes a reference, and the next is compared to it before updating it. Published manual runs read the latest scheduled snapshot (quick CPU uses the full CPU reference), but do not modify it. The PR comparison reads the CPU snapshot too. The snapshots contain configured observables and workload information, not particle arrays. This is a rolling reference, not an immutable scientific validation target: an intentional physics change will also become the next scheduled reference.

Manual runs of the scheduled reference, FLUKA, and Geant4 benchmark workflows upload their detailed artefacts but do not publish to permanent history by default. Enable the `publish_history` checkbox when dispatching a run to publish its compact record, compare to the scheduled physics snapshot, and redeploy the dashboard. Scheduled runs always publish on success (external engines must have measurements). The separate PR workflows publish automatically as described above. The accepted Everest reference is different from these rolling benchmark snapshots: it is committed to the source branch, under `tests/data/everest_physics_reference/`, and is only updated deliberately.

The website is regenerated from **all** permanent records. Scheduled profiles share one column only for the same host, UTC day, branch, tested commit, Xcoll release, and recorded common dependency versions (latest successful run per profile and configuration). Different configurations get separate columns, even on the same day. Engine-specific versions (FLUKA, LineBuilder, FlukaIO, Geant4, BDSIM) are retained per component rather than used to split compatible profiles. Published manual runs remain separate. Contributing revisions, timestamps, and versions are retained. History writes and deployments have separate concurrency queues. `xcoll-performance-publish.yml` is shared by the external engines; `xcoll-performance-dashboard.yml` deploys the existing history without running benchmarks. Pages must use **GitHub Actions** as its publishing source.

## GPU performance

The CUDA benchmarks are opt-in (`cd benchmarks; pytest -m gpu`). They cover
collimators, blocks, crystals, and InteractionRecord export. They are not run
by the scheduled CPU workflows.

To choose a particle count for a GPU job, run from the repository root:

```bash
python benchmarks/profile_gpu.py --workload block --particles 100 1000 10000 100000 1000000
python benchmarks/profile_gpu.py --workload lhc --particles 100 1000 10000 50000 --turns 3
```

Every run measures serial CPU, OpenMP (`omp_num_threads='auto'` by default),
and CUDA on the same initial particle distribution. Use `--cpu-threads N`
to change the OpenMP team size without removing the serial comparison. The
script reports the observed OpenMP team size, runtime maximum, and available
logical CPUs. Use `--device` to select CUDA, and `--beam`, `--plane`,
`--crystals`, and `--turns` for the LHC workload.

The script prints warmed medians and the first sampled count where CUDA stays
faster than each CPU configuration, both for tracking and for input copy +
tracking + loss-map analysis. Initial compilation, warm-up, and pencil
generation are excluded from the warmed timings. Line/element preparation and
excluded warm-up are reported separately. These are machine- and
workload-specific measurements, not universal thresholds. Repeated timing tables are printed together at the end, and the crossover plot is saved to a timestamped PNG in `benchmarks/results/`. Use `--plot PATH` to select the output location or `--no-plot` to disable plotting.

## Physics validation

The Everest reference-validation tests compare the candidate's CPU physics against an accepted reference committed in `tests/data/everest_physics_reference/`. This is distinct from the scalar physics observables stored with the performance benchmarks: it is a fixed validation target, not a rolling reference that is replaced on every successful run.

### Run the physics tests

From the repository root, with the candidate installed in editable mode and the test dependencies available:

```bash
python -m pip install -e '.[tests,benchmarks]'

# Rebuild the serial Xsuite kernels for this checkout (adjust the number of compiler workers).
python -m benchmarks.perf.prebuild_cpu \
    --skip-tpsa --skip-synrad --n-threads 12 --clear-existing

XCOLL_EVEREST_PHYSICS_FULL_MATRIX=1 \
XCOLL_EVEREST_PHYSICS_USE_PREBUILT=1 \
python -m pytest tests/test_everest_physics_reference.py -v --durations=0
```

No reference generation is needed. Without an explicit reference override, the test uses the committed dataset and checks the full matrix. The prebuild step compiles kernels against the current candidate; only enable `XCOLL_EVEREST_PHYSICS_USE_PREBUILT=1` after rebuilding for that exact checkout. Otherwise the tests force kernel compilation, to avoid using a prebuilt kernel containing an older implementation of Xcoll. Rebuild after changing kernel sources or relevant dependencies. Prebuilding excludes TPSA and synchrotron-radiation configurations, but auxiliary kernels not covered by the prebuild may still be compiled during the tests.

To run a subset, add a pytest expression, e.g. `-k "not lhc"`, `-k "lhc and b1"`, or `-k "block and tungsten"`. The full reference remains unchanged. Local uncommitted candidate changes are allowed; the imported Xcoll package must come from the checkout being tested.

The default matrix contains 40 physics cases:

- 24 standalone cases: CFC and tungsten blocks and collimators, and four silicon crystal regimes (channelling, volume reflection, amorphous transport, and miss), each at 26, 450, and 7000 GeV/c.
- 16 LHC cases: both beams, both planes, with and without crystals, at 450 and 7000 GeV/c. This is an energy scan of the same normalised Run 3 lattice, not a comparison of independently matched injection and flat-top optics.

The tests replay the stored initial coordinates and RNG states exactly, including on a different machine. They use the reference's recorded particle populations, momenta, tracking turns, and other workload parameters. The reference manifest records the source commit, Python, dependencies, and platform. Environment differences are reported but do not themselves invalidate a comparison; incompatible physical inputs, workload parameters, or LHC data do. The physics PR workflow installs the reference's pinned Python and dependencies, apart from the editable candidate Xcoll package.

### What is compared?

The comparison checks survival and particle-state fractions, survivor coordinate, energy, and kick distributions (with sample-size-aware KS screening limits), and channelled fractions where applicable. For the standalone elements we also compare per-primary interaction multiplicities, incidences, and the total longitudinal extent of continuous interactions. The paired statistics use differences between the same incident primaries, rather than treating repeated interactions as independent particles. For continuous processes, comparing the accumulated extent makes the check insensitive to how a trajectory is split into recorded segments. Deposited energy and, for the LHC cases, collimator and aperture loss distributions are checked as well.

The tolerances are screening criteria, not validated physics acceptance thresholds. A failed comparison requires investigation: do not replace the accepted reference simply to make a candidate pass. The standalone comparison script prints the individual measured differences and limits; the PR dashboard preserves these results from the workflow reports.

### Generate or update an accepted reference

Use [everest_compare_physics.py](everest_compare_physics.py), not the former `compare_everest_physics.py` script. A new permanent reference should be generated only from an explicitly validated revision, with the full matrix and the intended environment. The current reference uses 20,000 particles for most standalone cases and all LHC cases, ten LHC turns, and higher populations for selected rare-event workloads.

For example, on a clean checkout of the validated revision, rebuild the kernels as above and run:

```bash
ts=$(date "+%Y%m%d-%H%M%S")
python benchmarks/everest_compare_physics.py collect \
    --checkout "$PWD" --use-prebuilt \
    --output "/tmp/everest-physics-$ts" \
    --particles 20000 --lhc-particles 20000 --turns 10 \
    --events-per-particle 12 \
    --case-particles block-tungsten-26gev=500000 \
    --case-particles crystal-amorphous-26gev=100000 \
    --case-particles block-cfc-7000gev=100000
```

Validate the generated dataset before replacing the committed reference:

```bash
XCOLL_EVEREST_PHYSICS_REFERENCE="/tmp/everest-physics-$ts" \
XCOLL_EVEREST_PHYSICS_FULL_MATRIX=1 \
XCOLL_EVEREST_PHYSICS_USE_PREBUILT=1 \
python -m pytest tests/test_everest_physics_reference.py -v --durations=0
```

Inspect and review the physics results, including the comparison with the previously accepted dataset where the workload remains compatible. Then deliberately replace the entire `tests/data/everest_physics_reference/` directory, including `reference.json`, the case NPZ files, and `inputs/arrays.npz`, and commit the new reference. Keep ordinary generated results in temporary or ignored locations. The compact storage shares identical initial-particle chunks between cases and losslessly encodes the inputs; the accepted reference must be retained as a complete dataset.

### Compare other revisions or existing snapshots

The `reference` command resolves a branch, tag, or commit to a full SHA, creates and removes a temporary detached worktree, and leaves the current checkout unchanged. It deliberately forces compilation from that revision; for efficient reference generation on a validated checkout, prefer the explicit prebuild and `collect --use-prebuilt` recipe above.

The `collect` and `reference` commands accept `--input-reference` to replay exactly the initial particles and settings of an existing reference. For example, to compare the current editable checkout against the committed reference, after rebuilding its kernels:

```bash
python benchmarks/everest_compare_physics.py collect \
    --checkout "$PWD" --use-prebuilt \
    --input-reference tests/data/everest_physics_reference \
    --output /tmp/physics-candidate

python benchmarks/everest_compare_physics.py compare \
    tests/data/everest_physics_reference /tmp/physics-candidate
```

Choose a new output directory for each collection. `XCOLL_EVEREST_PHYSICS_REFERENCE` selects an alternative reference for pytest, and `XCOLL_EVEREST_PHYSICS_REF_SHA` can enforce its expected source SHA. External references may contain a subset of the matrix unless `XCOLL_EVEREST_PHYSICS_FULL_MATRIX=1` is set. See the module docstring and `--help` for case and energy selection, separate worktrees, and the lossless `compact` command for existing schema-4 or schema-5 snapshots.

### Physics validation in CI

Both automatic and manually dispatched physics PR runs load the accepted snapshots from a separately pinned main checkout and run the current validation harness against the pinned PR candidate. The candidate's Xcoll source and lattice data are used for tracking; the harness's latest snapshot decoder reads the reference. Each of the three hosted matrix jobs (standalone elements, LHC B1, and LHC B2) rebuilds its candidate serial kernels independently. A missing accepted reference is an error, not a reason to generate one in CI.

The tested base, head, merge, validation-tooling SHA, and accepted reference-source SHA are recorded separately. This allows a physics run to validate an existing historical-base PR merge without requiring a branch update or restarting the CPU performance benchmarks. The separate publisher saves the outcomes and deploys the dashboard automatically. See [PR_WORKFLOWS.md](PR_WORKFLOWS.md) for the manual dispatch commands and the publication details.

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
