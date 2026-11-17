"""Display-only benchmark presentation regressions."""

import json
from pathlib import Path

from perf import results


def test_compact_filename_and_unchanged_resolution(tmp_path, capsys):
    storage = tmp_path / 'results'
    platform = storage / 'Darwin-CPython-3.14-64bit'
    platform.mkdir(parents=True)
    sha = 'a' * 40
    filename = f'0027_{sha}_20261117_124500.json'
    path = platform / filename
    path.write_text(json.dumps({
        'datetime': '2026-11-17T12:45:00',
        'ci': {'source_branch': 'main', 'source_sha': sha},
        'benchmarks': [{'name': 'sample'}],
    }))

    assert results.main(['--results', str(storage), 'list']) == 0
    output = capsys.readouterr().out
    assert 'Darwin-CPython-3.14-64bit/0027_aaaaaaaa….json' in output
    assert filename not in output
    assert results.resolve('1', storage, results.entries(storage), {}) == path
    # Display abbreviations must not affect archived or tagged paths.
    archive = storage / 'archive' / 'snapshot' / platform.name
    archive.mkdir(parents=True)
    archived = archive / filename
    archived.write_text(path.read_text())
    assert results._display_filename(archived, storage).startswith('archive/snapshot/')


def test_nonstandard_names_and_physics_marker(tmp_path):
    storage = tmp_path / 'results'
    path = storage / 'MacOS' / '0002_sample.json'
    assert results._display_filename(path, storage) == 'MacOS/0002_sample.json'
    dashboard = (Path(__file__).resolve().parents[1] / 'perf' / 'templates'
                 / 'dashboard.html').read_text(encoding='utf-8')
    assert '#timeline svg path[fill="#b64030"]' in dashboard
    assert 'fill:#ffe075' in dashboard
    assert "['#ffe075','physics changed (gold corner)']" in dashboard
