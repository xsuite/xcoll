"""Git provenance regression tests (no Xsuite dependencies)."""

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from perf.provenance import git_state, check_state


@pytest.fixture
def repo(tmp_path):
    def git(*args):
        subprocess.check_call(['git', '-C', str(tmp_path), *args],
                              stdout=subprocess.DEVNULL)
    git('init', '-q')
    git('config', 'user.name', 'Benchmark Test')
    git('config', 'user.email', 'test@example.invalid')
    (tmp_path / 'code.py').write_text('x = 1\n')
    (tmp_path / '.gitignore').write_text('results/\n')
    git('add', '.')
    git('commit', '-qm', 'initial')
    return tmp_path, git


def test_unchanged_dirty_checkout(repo):
    root, _ = repo
    (root / 'code.py').write_text('x = 2\n')
    (root / 'new.py').write_text('y = 1\n')
    start = git_state(root)
    (root / 'results').mkdir()
    (root / 'results' / 'run.json').write_text('{}')
    check_state(start, git_state(root), start.commit)
    (root / 'new.py').write_text('y = 2\n')
    with pytest.raises(ValueError, match='provenance error'):
        check_state(start, git_state(root), start.commit)


def test_dirty_to_dirty_and_staging(repo):
    root, git = repo
    (root / 'code.py').write_text('x = 2\n')
    start = git_state(root)
    (root / 'code.py').write_text('x = 3\n')
    with pytest.raises(ValueError):
        check_state(start, git_state(root), start.commit)
    start = git_state(root)
    git('add', 'code.py')
    with pytest.raises(ValueError):
        check_state(start, git_state(root), start.commit)


def test_commit_and_branch_change(repo):
    root, git = repo
    start = git_state(root)
    git('checkout', '-qb', 'other')
    with pytest.raises(ValueError):
        check_state(start, git_state(root), start.commit)
    start = git_state(root)
    git('commit', '--allow-empty', '-qm', 'new commit')
    with pytest.raises(ValueError):
        check_state(start, git_state(root), git_state(root).commit)


def test_reported_commit_mismatch(repo):
    root, _ = repo
    start = git_state(root)
    with pytest.raises(ValueError):
        check_state(start, start, '0' * 40)