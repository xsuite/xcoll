"""Check that a benchmark's checkout is unchanged before publishing results."""

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import os
import subprocess


@dataclass(frozen=True)
class GitState:
    commit: str
    branch: str
    fingerprint: str


def git_state(root):
    root = Path(root)

    def git(*args):
        return subprocess.check_output(
            ['git', '-C', str(root), *args]
        )

    commit = git('rev-parse', 'HEAD').decode().strip()
    branch = git('rev-parse', '--abbrev-ref', 'HEAD').decode().strip()
    h = sha256()

    def add(data):
        h.update(len(data).to_bytes(8, 'big'))
        h.update(data)

    # Capture both the index and working tree,
    # including binary changes.
    add(git(
        'diff', '--no-ext-diff', '--no-textconv',
        '--binary', 'HEAD', '--'
    ))
    add(git(
        'diff', '--no-ext-diff', '--no-textconv',
        '--binary', '--cached', 'HEAD', '--'
    ))

    # New, non-ignored source files can also affect
    # a dirty-checkout run.
    untracked = git(
        'ls-files', '--others', '--exclude-standard', '-z'
    ).split(b'\0')

    for name in sorted(untracked):
        if not name:
            continue

        path = root / os.fsdecode(name)
        add(name)

        if path.is_symlink():
            add(b'link')
            add(os.fsencode(os.readlink(path)))
        elif path.is_file():
            add(b'file')
            add(str(path.stat().st_mode).encode())
            add(path.read_bytes())
        else:
            raise RuntimeError(
                f'Cannot fingerprint untracked path: {path}'
            )

    return GitState(commit, branch, h.hexdigest())


def check_state(start, end, reported_commit):
    if start != end or reported_commit != start.commit:
        raise ValueError(
            'Benchmark provenance error: checkout changed '
            'or Git metadata disagrees.\n'
            f'Start: {start.commit} ({start.branch})\n'
            f'End:   {end.commit} ({end.branch})\n'
            f'JSON commit: {reported_commit}\n'
            'Results have not been saved. '
            'Rerun on an unchanged checkout.'
        )
