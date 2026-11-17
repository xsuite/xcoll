"""Regression checks for local report comparisons and dynamic workflow names."""
import json

import pytest

from benchmarks.perf import pr_extra, report


def _measurement(path, timestamp, lost, *, particles=1_000_000):
    record = {
        'datetime': timestamp,
        'machine_info': {'system': 'Darwin', 'processor': 'test CPU'},
        'commit_info': {'branch': 'main', 'id': 'a' * 40},
        'benchmarks': [{
            'name': 'bench_collimator_cfc_hit_1M',
            'fullname': 'bench_everest.py::bench_collimator_cfc_hit_1M',
            'group': 'Everest CPU | collimator',
            'stats': {'median': .01, 'stddev': 0, 'rounds': 5},
            'extra_info': {
                'protocol': 2, 'particles': particles,
                'physics_case': 'collimator-CFC-hit',
                'lost_fraction': lost,
            },
        }],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record), encoding='utf-8')
    return path


def test_local_report_numbers_tags_and_physics_comparison(tmp_path):
    results = tmp_path / 'results'
    older = _measurement(results / 'a' / '0001_baseline.json',
                         '2026-11-12T10:00:00', .675234)
    newer = _measurement(results / 'a' / '0002_candidate.json',
                         '2026-11-13T10:00:00', 0.)
    (results / '.labels').write_text(json.dumps({'baseline': str(older.relative_to(results))}))
    original = (older.read_bytes(), newer.read_bytes())

    # Run 1 is newest. The local report compares chronologically, even when
    # the user enters newest first, and marks the *newer* run as CHANGED.
    output = tmp_path / 'dashboard.html'
    assert report.main(['--results', str(results), '1', 'baseline', '-o', str(output)]) == 0
    html = output.read_text(encoding='utf-8')
    assert '"status":"CHANGED"' in html
    assert 'lost_fraction' in html
    assert 'One or more physical metrics outside tolerance' in html
    assert (older.read_bytes(), newer.read_bytes()) == original


def test_local_report_one_file_and_directory_remain_supported(tmp_path):
    results = tmp_path / 'results'
    case = _measurement(results / '0001_result.json', '2026-11-12T10:00:00', .675234)
    output = tmp_path / 'single.html'
    assert report.main(['--results', str(results), '1', '-o', str(output)]) == 0
    assert '"status":"CHANGED"' not in output.read_text(encoding='utf-8')
    assert report.main([str(results), '-o', str(output)]) == 0
    assert output.is_file()
    # Existing path support must not require a run number.
    assert report.main([str(case), '-o', str(output)]) == 0


def test_local_report_incompatible_workloads_are_not_changes(tmp_path):
    old = _measurement(tmp_path / 'older.json', '2026-11-12T10:00:00', .675234)
    new = _measurement(tmp_path / 'newer.json', '2026-11-13T10:00:00', .0,
                       particles=100_000)
    runs = report.local_comparison_runs([new, old], report.load_config())
    assert runs[-1]['entries'][0]['physics']['status'] == 'INCOMPATIBLE'


@pytest.mark.parametrize('engine', ['fluka', 'geant4'])
def test_publisher_accepts_dynamic_external_run_names(tmp_path, monkeypatch, engine):
    run = {'name': f'{engine} benchmark for PR #241',
           'path': '.github/workflows/xcoll-performance-external-pr.yml',
           'repository': {'full_name': pr_extra.REPO},
           'id': 12345, 'run_attempt': 1}
    assert pr_extra.source_workflow(run) == 'external'
    monkeypatch.setattr(pr_extra, 'event', lambda: {'workflow_run': run})
    monkeypatch.setattr(pr_extra, 'download', lambda *args, **kwargs: False)
    monkeypatch.setenv('RUNNER_TEMP', str(tmp_path))
    # Missing artifact means skipped publication, not an unexpected workflow.
    pr_extra.publish()


def test_publisher_rejects_wrong_source_even_with_matching_display_name():
    run = {'name': 'Xcoll external performance (PR)',
           'path': '.github/workflows/untrusted.yml',
           'repository': {'full_name': pr_extra.REPO}}
    with pytest.raises(ValueError, match='Unexpected source workflow'):
        pr_extra.source_workflow(run)
    run['path'] = '.github/workflows/xcoll-performance-external-pr.yml'
    run['repository']['full_name'] = 'another/repository'
    with pytest.raises(ValueError, match='Unexpected source repository'):
        pr_extra.source_workflow(run)
