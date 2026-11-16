#!/usr/bin/env bash
# Trusted main orchestration; only administrator-reviewed PR code may run here.
set -euo pipefail
case "$ENGINE" in fluka|geant4) ;; *) exit 2 ;; esac
case "$COMPILE_WORKERS" in 0|2) ;; *) exit 2 ;; esac
work="$RUNNER_TEMP/xcoll-external-pr-${GITHUB_RUN_ID}-${GITHUB_RUN_ATTEMPT}"
mkdir -p "$work"
base="$GITHUB_WORKSPACE/baseline"
candidate="$GITHUB_WORKSPACE/candidate"
control="$GITHUB_WORKSPACE/control"
ctx="$RUNNER_TEMP/xcoll-external-pr-context-${GITHUB_RUN_ID}-${GITHUB_RUN_ATTEMPT}/pr-context.json"
mamba="$MINIFORGE_ROOT/bin/mamba"
export BENCHMARK_TESTBED="cern-openstack-xcoll-${ENGINE}-pr-paired"
# Prevent state reuse; the shared CERN concurrency group excludes daily jobs.
set +u
source "$MINIFORGE_ROOT/etc/profile.d/conda.sh"
set -u
activate() {
    set +u
    conda activate "$1"
    set -u
    test "$CONDA_PREFIX" = "$1"
}
if [ "$ENGINE" = fluka ]; then
    packages=(python=3.12 pip make cmake meson ninja gcc_linux-64 gxx_linux-64 gfortran_linux-64)
    extra=(cython)
    command -v gfortran >/dev/null
    : "${XCOLL_FLUKA_ROOT:?Missing FLUKA root}" "${XCOLL_FLUKA_SERVER:?Missing FLUKA server}" \
      "${XCOLL_FLUKAIO_LIB:?Missing FlukaIO}" "${XCOLL_LINEBUILDER_ROOT:?Missing LineBuilder}"
else
    packages=(python=3.14 pip bdsim-g4 make cmake 'hdf5>=2,<3=nompi_*' 'expat=2.8.1' cxx-compiler)
    extra=(cython pybind11 rpyc)
fi
"$mamba" create --yes --prefix "$work/base-env" --override-channels --channel conda-forge "${packages[@]}"
activate "$work/base-env"
"$mamba" list --prefix "$CONDA_PREFIX" --explicit > "$work/conda-spec.txt"
python -m pip install --upgrade pip
python -m pip install --upgrade xsuite
(cd "$base"; python -m pip install -e '.[benchmarks]' "${extra[@]}")
python -m pip freeze --exclude xcoll > "$work/pip-lock.txt"
# Record and reuse the exact package solve, not a second independent solve.
run_one() {
    local tree="$1" label="$2" env="$3"
    activate "$env"
    export XCOLL_PATH="$work/$label-state"
    cd "$tree"
    rm -rf benchmarks/results
    python -m benchmarks.perf.prebuild_cpu --skip-tpsa --skip-synrad --n-threads "$COMPILE_WORKERS" --clear-existing
    python -m benchmarks.perf.setup_external "$ENGINE"
    python -m pip freeze > "$work/$label-dependencies.txt"
    python -m benchmarks.perf.dependencies >> "$work/$label-dependencies.txt"
    (cd benchmarks; python -m pytest -m "$ENGINE" --basetemp="$work/$label-pytest")
    # pr_compare also rejects empty or multiple benchmark JSON files.
    test -n "$(find benchmarks/results -name '*.json' -print -quit)"
}
run_one "$base" base "$work/base-env"
"$mamba" create --yes --prefix "$work/candidate-env" --file "$work/conda-spec.txt"
activate "$work/candidate-env"
python -m pip install -r "$work/pip-lock.txt"
(cd "$candidate"; python -m pip install -c "$work/pip-lock.txt" -e '.[benchmarks]' "${extra[@]}")
run_one "$candidate" candidate "$work/candidate-env"
# Compare from trusted orchestration using the base environment. Never recompute
# in the history publisher, and never publish through a write token on CERN.
activate "$work/base-env"
cd "$control"
branch=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["branch"])' "$ctx")
python -m benchmarks.perf.pr_compare \
    --profile "$ENGINE" --base "$base/benchmarks/results" --candidate "$candidate/benchmarks/results" \
    --base-sha "$BASE_SHA" --candidate-sha "$CANDIDATE_SHA" --head-sha "$HEAD_SHA" \
    --pr-number "$PR_NUMBER" --branch "$branch" --run-id "$GITHUB_RUN_ID" \
    --attempt "$GITHUB_RUN_ATTEMPT" --json-output "$work/pr-comparison.json" \
    --output "$work/pr-comparison.md"
cat "$work/pr-comparison.md" >> "$GITHUB_STEP_SUMMARY"
