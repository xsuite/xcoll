"""Regression tests for portable paired-run pip locks."""
import json

import pytest

from benchmarks.perf import pip_lock


CONDA = [
    {'name': 'cffi', 'version': '2.0.0', 'channel': 'conda-forge'},
    {'name': 'numpy', 'version': '2.5.5', 'channel': 'conda-forge'},
    {'name': 'PyYAML', 'version': '6.0.3', 'channel': 'conda-forge'},
]
VERSIONS = {'cffi': '2.0.0', 'numpy': '2.5.5', 'pyyaml': '6.0.3',
            'xsuite': '0.68.0', 'xtrack': '0.118.7'}


def installed_version(name):
    return VERSIONS[pip_lock.normalise(name)]


def test_conda_local_build_paths_not_replayed():
    frozen = '''cffi @ file:///home/conda/feedstock_root/build_artifacts/cffi/work
numpy==2.5.5
PyYAML @ file:///home/conda/feedstock_root/build_artifacts/pyyaml/work
xsuite==0.68.0
xtrack==0.118.7
'''
    assert pip_lock.portable_lock(CONDA, frozen, installed_version=installed_version) == (
        'xsuite==0.68.0\nxtrack==0.118.7\n')


def test_upgraded_conda_distribution_kept_at_pip_version():
    upgraded = dict(VERSIONS, numpy='2.6.0')
    lock = pip_lock.portable_lock(CONDA, 'numpy==2.6.0\nxsuite==0.68.0\n',
                                 installed_version=lambda name: upgraded[pip_lock.normalise(name)])
    assert lock == 'numpy==2.6.0\nxsuite==0.68.0\n'


@pytest.mark.parametrize('frozen', [
    'unknown @ file:///home/conda/feedstock_root/work\n',
    'cffi @ file:///home/conda/feedstock_root/work\nxsuite==0.68.0\n',
    '-e /tmp/xcoll-helper\n',
    'xsuite==0.68.0\nxsuite==0.68.0\n',
])
def test_nonportable_or_duplicate_pip_sources_rejected(frozen):
    versions = dict(VERSIONS, cffi='3.0.0')
    with pytest.raises(ValueError, match='Non-portable|Unsupported|Duplicate'):
        pip_lock.portable_lock(CONDA, frozen,
                               installed_version=lambda name: versions.get(pip_lock.normalise(name), '1'))


def test_empty_conda_inventory_rejected():
    with pytest.raises(ValueError, match='Empty Conda'):
        pip_lock.portable_lock([], 'xsuite==0.68.0\n', installed_version=installed_version)


def test_no_pip_additions_rejected():
    with pytest.raises(ValueError, match='Empty pip-only lock'):
        pip_lock.portable_lock(CONDA, 'cffi==2.0.0\n', installed_version=installed_version)


def test_cli_writes_pip_only_lock(tmp_path, monkeypatch):
    inventory = tmp_path / 'conda.json'
    inventory.write_text(json.dumps(CONDA))
    output = tmp_path / 'pip-lock.txt'
    constraints = tmp_path / 'pip-constraints.txt'
    def command(args, **kw):
        if 'freeze' in args:
            return 'cffi @ file:///home/conda/work\nxsuite==0.68.0\n'
        return json.dumps([{'name':'cffi','version':'2.0.0'},
                           {'name':'xsuite','version':'0.68.0'},
                           {'name':'xcoll','version':'0.13.1'}])
    monkeypatch.setattr(pip_lock.subprocess, 'check_output', command)
    monkeypatch.setattr(pip_lock.metadata, 'version', lambda name: VERSIONS[pip_lock.normalise(name)])
    pip_lock.main(['--conda-inventory', str(inventory), '--output', str(output),
                   '--constraints', str(constraints)])
    assert output.read_text() == 'xsuite==0.68.0\n'
    assert constraints.read_text() == 'cffi==2.0.0\nxsuite==0.68.0\n'


def test_constraints_preserve_conda_versions_but_skip_xcoll():
    inv = [{'name': 'cffi', 'version': '2.0.0'},
           {'name': 'xcoll', 'version': '0.13.1'},
           {'name': 'xsuite', 'version': '0.68.0'}]
    assert pip_lock.portable_constraints(inv) == 'cffi==2.0.0\nxsuite==0.68.0\n'
    with pytest.raises(ValueError, match='Malformed'):
        pip_lock.portable_constraints([{'name': 'bad', 'version': '1;rm -rf /'}])
