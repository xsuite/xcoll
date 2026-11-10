"""List, label, and compare pytest-benchmark results stored locally.

    python -m perf.results list
    python -m perf.results tag 1 release-baseline
    python -m perf.results compare release-baseline 2
    python -m perf.results untag release-baseline

Numbers refer to the most recent entries in `list` (1 = newest).
Labels are local aliases stored in results/.labels, not changes to JSON files.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .paths import results_dir

DEFAULT_RESULTS = results_dir()


def _labels_path(results):
    return results / '.labels'


def _load_labels(results):
    p = _labels_path(results)
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else {}


def _save_labels(results, labels):
    results.mkdir(parents=True, exist_ok=True)
    p = _labels_path(results)
    tmp = p.with_suffix('.tmp')
    tmp.write_text(json.dumps(labels, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    tmp.replace(p)


def _metadata(path):
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
    except (ValueError, OSError):
        return None
    if not isinstance(data, dict) or not isinstance(data.get('benchmarks'), list):
        return None
    ci, commit = data.get('ci') or {}, data.get('commit_info') or {}
    return {'path': path,
            'when': data.get('datetime') or 'unknown',
            'branch': ci.get('source_branch') or commit.get('branch') or 'unknown',
            'sha': ci.get('source_sha') or commit.get('id') or 'unknown',
            'n': len(data['benchmarks'])}


def entries(results):
    if not results.exists():
        return []
    out = []
    for p in results.rglob('*.json'):
        if not p.is_file():
            continue
        m = _metadata(p)
        if m:
            out.append(m)
    return sorted(out, key=lambda x: (x['when'], x['path'].name), reverse=True)


def resolve(selection, results, rows, labels):
    if selection.isdecimal():
        i = int(selection) - 1
        if i < 0 or i >= len(rows):
            raise ValueError(f'Run number {selection} is outside 1..{len(rows)}')
        return rows[i]['path']
    if selection in labels:
        p = results / labels[selection]
        if not p.is_file():
            raise ValueError(f'Tag {selection} points to missing result {p}')
        return p
    p = Path(selection)
    if not p.exists():
        p = results / selection
    if not p.is_file():
        raise ValueError(f'No result or tag called {selection!r}')
    return p


def archive_results(results, threshold=9000, force=False):
    """Rotate active autosaves before pytest-benchmark's 4-digit counter limit.

    Active storage is emptied, resetting its counter to 0001. All files remain
    available in a dated archive subdirectory; local tags are updated.
    """
    results = Path(results)
    labels = _load_labels(results)
    moved = []
    if not results.exists():
        return moved
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    for directory in [results, *[p for p in results.iterdir() if p.is_dir() and p.name != 'archive']]:
        files = [p for p in directory.glob('*.json')
                 if re.match(r'^\d{4,}_.+\.json$', p.name)]
        max_counter = max((int(p.name.split('_', 1)[0]) for p in files), default=0)
        if not force and len(files) < threshold and max_counter < threshold:
            continue
        for old in files:
            relative = old.relative_to(results)
            dest = results / 'archive' / stamp / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                raise FileExistsError(f'Refusing to overwrite archive {dest}')
            shutil.move(str(old), str(dest))
            old_rel, new_rel = str(relative), str(dest.relative_to(results))
            for tag, path in tuple(labels.items()):
                if path == old_rel:
                    labels[tag] = new_rel
            moved.append(dest)
    if moved:
        _save_labels(results, labels)
    return moved


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, default=DEFAULT_RESULTS)
    sub = p.add_subparsers(dest='command', required=True)
    sub.add_parser('list', help='Display saved runs (newest first)')
    archive = sub.add_parser('archive', help='Move old autosaves into archive/ without deleting them')
    archive.add_argument('--threshold', type=int, default=9000)
    archive.add_argument('--force', action='store_true', help='Archive the active autosaves now')
    tag = sub.add_parser('tag', help='Give a run a local readable alias')
    tag.add_argument('run', help='Row number, tag, or JSON file')
    tag.add_argument('name')
    untag = sub.add_parser('untag', help='Remove an alias without deleting the run')
    untag.add_argument('name')
    cmp = sub.add_parser('compare', help='Compare timings and configured observables')
    cmp.add_argument('baseline')
    cmp.add_argument('candidate')
    cmp.add_argument('--no-timings', action='store_true')
    cmp.add_argument('--allow-unverified', action='store_true')
    args = p.parse_args(argv)
    results = args.results.resolve()
    rows, labels = entries(results), _load_labels(results)
    if args.command == 'list':
        reverse_labels = {}
        for name, rel in labels.items():
            reverse_labels.setdefault(rel, []).append(name)
        print(f"{'#':>4}  {'Timestamp':19}  {'Branch':22}  {'Commit':10}  {'Cases':>5}  Tag                     Filename")
        print('-' * 135)
        for i, m in enumerate(rows, 1):
            rel = str(m['path'].relative_to(results))
            tags = ','.join(sorted(reverse_labels.get(rel, [])))
            print(f"{i:>4}  {m['when'][:19]:19}  {m['branch'][:22]:22}  {m['sha'][:10]:10}  {m['n']:>5}  {tags:22}  {rel}")
        print(f'({len(rows)} saved run(s) in {results})')
        return 0
    if args.command == 'archive':
        if args.threshold < 1:
            p.error('--threshold must be positive')
        moved = archive_results(results, threshold=args.threshold, force=args.force)
        print(f'Archived {len(moved)} result(s) under {results / "archive"}')
        return 0
    if args.command == 'tag':
        if not re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_.-]{0,63}', args.name):
            p.error('Tag must contain 1–64 letters, numbers, _, -, or ., starting alphanumeric')
        if args.name in labels:
            p.error(f'Tag already exists: {args.name}. Untag it first.')
        file = resolve(args.run, results, rows, labels).resolve()
        try:
            rel = str(file.relative_to(results))
        except ValueError:
            p.error('Only local files under the results directory can be tagged')
        labels[args.name] = rel
        _save_labels(results, labels)
        print(f'{args.name} -> {rel}')
        return 0
    if args.command == 'untag':
        if args.name not in labels:
            p.error(f'No such tag: {args.name}')
        del labels[args.name]
        _save_labels(results, labels)
        return 0
    if args.command == 'compare':
        from . import compare  # native pytest-benchmark CLI remains unchanged
        params = [str(resolve(args.baseline, results, rows, labels)),
                  str(resolve(args.candidate, results, rows, labels))]
        if args.no_timings:
            params.append('--no-timings')
        if args.allow_unverified:
            params.append('--allow-unverified')
        return compare.main(params)
    return 2


if __name__ == '__main__':
    raise SystemExit(main())
