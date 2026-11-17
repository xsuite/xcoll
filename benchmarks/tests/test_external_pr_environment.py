"""Exercise the paired external benchmark setup without Conda or an engine."""
import os
from pathlib import Path
import subprocess

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'perf' / 'run_external_pr.sh'
LOCK = ('# platform: linux-64\n@EXPLICIT\n'
        'https://conda.anaconda.org/conda-forge/linux-64/python-3.14-0.conda'
        '#0123456789abcdef0123456789abcdef\n')


@pytest.fixture
def fake_runner(tmp_path):
    root = tmp_path / 'miniforge'
    bin_dir = root / 'bin'
    bin_dir.mkdir(parents=True)
    profile = root / 'etc/profile.d/conda.sh'
    profile.parent.mkdir(parents=True)
    profile.write_text('conda() { test "$1" = activate; export CONDA_PREFIX="$2"; }\n')
    # Record actual commands and simulate only the external side effects of
    # the script. In particular, the real Bash lock validation is exercised.
    stub = r'''#!/usr/bin/env bash
set -euo pipefail
name=${0##*/}
{ printf '%s\t%s' "$name" "${CONDA_PREFIX:-}"; printf '\t%s' "$@"; printf '\n'; } >> "$CALLS"
case "$name" in
    conda)
        test "$1" = list
        test "$2" = --prefix
        case "$4" in
            --explicit)
                test "$5" = --md5
                printf '%s' "$SPEC"
                ;;
            --json)
                printf '%s\n' '[{"name":"cffi","version":"2.0.0","channel":"conda-forge"}]'
                ;;
            *) exit 2 ;;
        esac
        ;;
    mamba)
        test "$1" = create
        test "$3" = --prefix
        if [ "$#" -eq 6 ] && [ "$5" = --file ]; then
            test "${FAIL_CREATE:-0}" != 1
            test "$(cat "$6")" = "$(printf '%s' "$SPEC")"
        fi
        mkdir -p "$4"
        ;;
    python)
        if [[ "$1" = */benchmarks/perf/pip_lock.py ]]; then
            shift
            while (( $# )); do
                case "$1" in
                    --output)
                        printf 'xsuite==0.68.0\n' > "$2"
                        shift 2
                        ;;
                    --constraints)
                        printf 'cffi==2.0.0\nxsuite==0.68.0\n' > "$2"
                        shift 2
                        ;;
                    *) shift ;;
                esac
            done
        elif [ "$1" = -m ] && [ "$2" = pytest ]; then
            mkdir -p results
            printf '{}\n' > results/benchmark.json
        elif [ "$1" = -m ] && [ "$2" = benchmarks.perf.pr_compare ]; then
            while [ "$1" != --output ]; do shift; done
            printf 'comparison\n' > "$2"
        elif [ "$1" = -c ]; then
            printf 'feature\n'
        fi
        ;;
esac
'''
    for name in ('conda', 'mamba', 'python', 'gfortran'):
        path = bin_dir / name
        path.write_text(stub + '\n')
        path.chmod(0o755)
    workspace = tmp_path / 'workspace'
    for checkout in ('baseline', 'candidate', 'control'):
        (workspace / checkout / 'benchmarks').mkdir(parents=True)
    temp = tmp_path / 'runner-temp'
    temp.mkdir()
    calls = tmp_path / 'calls.jsonl'
    env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ['PATH'],
               MINIFORGE_ROOT=str(root), GITHUB_WORKSPACE=str(workspace),
               RUNNER_TEMP=str(temp), GITHUB_RUN_ID='42', GITHUB_RUN_ATTEMPT='1',
               GITHUB_STEP_SUMMARY=str(tmp_path / 'summary'), COMPILE_WORKERS='2',
               ENGINE='geant4', BASE_SHA='a'*40, HEAD_SHA='b'*40,
               CANDIDATE_SHA='c'*40, PR_NUMBER='241', CALLS=str(calls), SPEC=LOCK,
               XCOLL_FLUKA_ROOT='/fake/root', XCOLL_FLUKA_SERVER='/fake/server',
               XCOLL_FLUKAIO_LIB='/fake/lib', XCOLL_LINEBUILDER_ROOT='/fake/linebuilder')

    def run(**overrides):
        result = subprocess.run(['bash', str(SCRIPT)], env=dict(env, **overrides),
                                capture_output=True, text=True, cwd=workspace)
        trace = []
        for line in calls.read_text().splitlines():
            name, prefix, *args = line.split('\t')
            trace.append((name, args, prefix))
        return result, trace

    return run


@pytest.mark.parametrize('engine', ['fluka', 'geant4'])
def test_exact_conda_lock_created_before_benchmark(fake_runner, engine):
    result, trace = fake_runner(ENGINE=engine)
    assert result.returncode == 0, result.stderr
    assert [x[0] for x in trace[:4]] == ['mamba', 'conda', 'conda', 'mamba']
    assert trace[2][1][3] == '--json'
    assert '--file' in trace[3][1]
    assert trace[1][1][2].endswith('/base-env')
    assert trace[3][1][3].endswith('/candidate-env')
    lock_command = next(i for i, x in enumerate(trace) if x[0] == 'python'
                        and x[1] and x[1][0].endswith('/benchmarks/perf/pip_lock.py'))
    tracking_indices = [i for i, x in enumerate(trace) if x[0] == 'python'
                        and x[1][:2] == ['-m', 'pytest']]
    assert lock_command < tracking_indices[0]
    tracking = [trace[i] for i in tracking_indices]
    assert len(tracking) == 2
    assert tracking[0][2].endswith('/base-env')
    assert tracking[1][2].endswith('/candidate-env')
    assert all(x[1][3] == engine for x in tracking)


@pytest.mark.parametrize('spec', [
    '# platform: linux-64\n',
    '@EXPLICIT\n',
    '/home/xcollrunner/base-env\n' + LOCK,
    LOCK + '/home/xcollrunner/base-env\n',
])
def test_invalid_conda_lock_fails_before_compilation(fake_runner, spec):
    result, trace = fake_runner(SPEC=spec)
    assert result.returncode != 0
    assert 'Invalid Conda explicit package lock' in result.stderr
    assert [x[0] for x in trace] == ['mamba', 'conda', 'conda']


def test_candidate_create_failure_prevents_baseline_run(fake_runner):
    result, trace = fake_runner(FAIL_CREATE='1')
    assert result.returncode != 0
    assert [x[0] for x in trace] == ['mamba', 'conda', 'conda', 'mamba']
