"""Tests for benchmark machine labels and Xcoll release metadata."""
from __future__ import annotations

import json
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))
from perf import machine, report


def _saved_run(tmp_path, *, name='BlackMetal3', version='0.13.0', profile='full',
               when='2026-11-10T04:00:00', triggered='schedule'):
    path = tmp_path / f'{profile}.json'
    data = {
        'datetime': when,
        'machine_name': name,
        'machine_info': {
            'system': 'Darwin', 'processor': 'arm',
            'node': 'blackmetal3.some-random-bonjour-network.local',
        },
        'xsuite_versions': {'xcoll': version},
        'testbed': 'cern-openstack-xcoll-cpu-v1',
        'ci': {
            'source_branch': 'main', 'source_sha': 'a' * 40,
            'profile': profile, 'trigger': triggered,
            'physical_host': 'cern-openstack-xcoll-vm-v1',
        },
        'benchmarks': [
            {'name': f'bench_{profile}', 'stats': {'median': 0.01},
             'extra_info': {'protocol': 2}}
        ],
    }
    path.write_text(json.dumps(data))
    return path


def test_macos_uses_computer_name_not_network_name(monkeypatch):
    monkeypatch.delenv('BENCHMARK_MACHINE_NAME', raising=False)
    monkeypatch.setattr(machine.platform, 'system', lambda: 'Darwin')
    monkeypatch.setattr(machine, '_command', lambda *parts: 'BlackMetal3'
                        if parts == ('scutil', '--get', 'ComputerName') else None)
    assert machine.machine_name() == 'BlackMetal3'


def test_machine_override_is_honoured(monkeypatch):
    monkeypatch.setenv('BENCHMARK_MACHINE_NAME', 'BlackMetal3')
    monkeypatch.setattr(machine.platform, 'system', lambda: 'Linux')
    assert machine.machine_name() == 'BlackMetal3'


def test_machine_fallback_does_not_use_bonjour_node(monkeypatch):
    monkeypatch.delenv('BENCHMARK_MACHINE_NAME', raising=False)
    monkeypatch.setattr(machine.platform, 'system', lambda: 'Darwin')
    monkeypatch.setattr(machine.platform, 'machine', lambda: 'arm64')
    monkeypatch.setattr(machine, '_command', lambda *parts: None)
    assert machine.machine_name() == 'Darwin | arm64'


def test_report_records_release_and_computer_name(tmp_path):
    p = _saved_run(tmp_path)
    (run,) = report.read_runs([p])
    assert run['machine'] == 'BlackMetal3'
    assert run['release'] == '0.13.0'
    assert run['components'][0]['release'] == '0.13.0'


def test_online_daily_release_mismatch_and_machine_name(tmp_path):
    files = [
        _saved_run(tmp_path, profile='full', version='0.13.0'),
        _saved_run(tmp_path, profile='fluka', version='0.13.1',
                   when='2026-11-10T06:00:00'),
    ]
    daily = report.read_runs(files, online=True)
    assert len(daily) == len({r['id'] for r in daily}) == 2
    assert all(r['machine'] == 'BlackMetal3' for r in daily)
    assert {r['release'] for r in daily} == {'0.13.0', '0.13.1'}
    assert all(len(r['components']) == 1 for r in daily)


def test_dashboard_header_has_four_labels():
    from pathlib import Path
    template = (Path(__file__).resolve().parents[1] / 'perf' /
                'templates' / 'dashboard.html').read_text()
    assert '`Xcoll ${sec.release||' in template
    assert 'y:72' in template
    assert 'r.machine||r.physical_host' in template
