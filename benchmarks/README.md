# Benchmarks

Xcoll benchmarks use [pytest-benchmark](https://pytest-benchmark.readthedocs.io/). They run serially, store individual timing samples, and record the Git revision and physical diagnostics. The benchmark functions are in `bench_*.py`; the supporting tools are in `perf/`.

## Running benchmarks

From the Xcoll repository root:

```bash
cd benchmarks
pytest                         # All benchmarks
pytest -m quick                # Quick profile
pytest -m full                 # Full profile only
pytest -m fluka                # Requires compiled FLUKA interface
pytest -m geant4               # Requires compiled Geant4/BDSIM interface
pytest -k crystal              # Benchmark names containing "crystal"
pytest --collect-only -q       # List available benchmarks
```

Each invocation saves a timestamped JSON file in `results/`. The directory is ignored by Git. Passing `-n` or `--tx` is disallowed; concurrent workers would invalidate timing measurements. Kernel compilation and preparation should occur outside the timed section.

## Saved results

From `benchmarks/`:

```bash
python -m perf.results list
python -m perf.results tag 1 release-013
python -m perf.results tag 2 gpu-branch
python -m perf.results compare release-013 gpu-branch
python -m perf.results untag gpu-branch
python -m perf.results archive
```

`list` prints the newest runs first, including timestamps, branches, and commits. Numbers refer to rows in the listing, with `1` denoting the newest. Tags are local aliases in `results/.labels`; JSON files are not renamed or modified. The `archive` command moves older files into `results/archive/` if the pytest-benchmark numbering approaches its discovery limit. Use `archive --force` to rotate the active results immediately. The archived results and tags remain accessible.

The `compare` command runs pytest-benchmark's standard timing comparison followed by checks of the saved physical diagnostics. These checks use `metrics.toml` and do not replace physics regression tests. For timing comparisons alone, use `pytest-benchmark compare` directly.

## Dashboard

```bash
python -m perf.report                    # All local results
python -m perf.report results -o results/dashboard.html
python -m perf.report baseline.json candidate.json -o comparison.html
```

Open the generated HTML file in a browser. It is self-contained and does not require a server or external JavaScript libraries. The default view shows `main` measurements as a chronological heatmap. Select a benchmark to display its runtime history, or use the comparison selectors to compare two stored revisions. Filters are available for testbed, branch, date, and benchmark name. Revisions from any branch or testbed can be selected; only identical workloads measured on the same testbed can be compared directly.

## Continuous integration

The following workflows are defined in `.github/workflows/`:

- `xcoll-performance.yml` runs the full profile daily at 03:17 UTC on `main`, using the dedicated CERN OpenStack runner. A manual run on `main` also uses that runner. A manually selected non-`main` revision runs on GitHub-hosted hardware instead. Successful measurements are committed to `benchmark-results` and used to regenerate GitHub Pages.
- `xcoll-performance-pr.yml` runs the quick profile for non-draft PRs on disposable GitHub-hosted runners. Its results are downloadable artifacts, but are not added to permanent history automatically.
- `xcoll-performance-fluka.yml` and `xcoll-performance-geant4.yml` run the corresponding external engine benchmarks manually on the restricted CERN runner. They rebuild their Xcoll interfaces in a disposable Conda environment, record dependencies, and publish results to the same history and dashboard. Schedules can be enabled after the installations are validated.
- `xcoll-performance-publish.yml` is reusable publication tooling shared by the two external workflows; it does not independently trigger a run.

The dedicated machine is selected by the label `xcoll-bench-cpu` and must have Miniforge installed at `/home/xcollrunner/miniforge3`. See [RUNNERS.md](RUNNERS.md) for setup and security requirements. GitHub Pages must use **GitHub Actions** as its publishing source, and the publishing job needs permission to write to `benchmark-results`.

The history branch stores immutable JSON records. The generated dashboard is published only after a successful measurement. GitHub-hosted runs should not be treated as a stable CPU reference across different physical hosts.

## Adding a benchmark

Add a standalone `bench_*` function to a `bench_*.py` file. Apply `@pytest.mark.quick`, `@pytest.mark.full`, `@pytest.mark.fluka`, or `@pytest.mark.geant4` and, where useful, `@pytest.mark.benchmark(group='...')`. Reset mutable particles, RNG state, and element counters between timed rounds. Record workload parameters and physical diagnostics in `benchmark.extra_info`.

The project-specific comparison contract lives in `metrics.toml`. Add the workload fields and observables there, not to the Python comparison or dashboard code. The modules in `perf/` contain no Xcoll-specific metric names; another Xsuite package can reuse them with its own configuration.

## Test the tooling

From the repository root:

```bash
python -m pytest -c /dev/null benchmarks/tests/test_tools.py -q
```

The tooling tests do not require an Xsuite physics installation.

## External engines and additional CPU measurements

`bench_import.py` measures `import xcoll` from a fresh Python interpreter. Its
standard pytest-benchmark median includes process start-up; the child-reported
**import-only** median is recorded separately as `child_import_median_s`. The
import is not a cold-disk-cache measurement.

`bench_materials.py`, `bench_rf_sweep.py`, `bench_blowup.py`, and
`bench_emittance_monitor.py` cover material construction/conversion, RF-sweep
preparation/tracking, random kicks, and monitor accumulation/reduction.
The external `bench_fluka.py` and `bench_geant4.py` files measure collimator
tracking with the engine already running and with fresh input particles for
every measured round. They deliberately exclude engine start-up, dependency
installation, and interface compilation.
