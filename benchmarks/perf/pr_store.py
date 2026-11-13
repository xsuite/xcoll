
"""Validate and retain compact PR benchmark comparisons from untrusted artifacts.

Only run this module from a trusted workflow checked out from main. The input
is treated exclusively as data, never as code, shell commands, or file paths.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

MAX_BYTES = 3_000_000
MAX_ROWS = 2000
SHA = re.compile(r'^[0-9a-f]{40}$')
ALLOWED_STATUSES = {'OK', 'CHANGED', 'INCOMPATIBLE', 'UNVERIFIED', 'NOT_CHECKED', 'BASELINE', 'MISSING'}


def _text(value, maximum=512):
    if not isinstance(value, str) or len(value) > maximum or any(ord(c) < 32 for c in value):
        raise ValueError('Invalid text field in PR comparison')
    return value


def _number(value):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise ValueError('Non-finite or invalid timing value')
    return value


def _physics(value):
    if value is None:
        return None
    if not isinstance(value, dict) or value.get('status') not in ALLOWED_STATUSES:
        raise ValueError('Invalid physics status')
    summary = _text(value.get('summary', ''), 500)
    details = value.get('details', [])
    if not isinstance(details, list) or len(details) > 50:
        raise ValueError('Invalid physics details')
    return {'status': value['status'], 'summary': summary,
            'details': [_text(x, 500) for x in details]}


def validate(record, *, run_id, attempt):
    if not isinstance(record, dict) or record.get('schema') != 1:
        raise ValueError('Unsupported PR comparison schema')
    pr = record.get('pr')
    if type(pr) is not int or not 0 < pr < 1_000_000:
        raise ValueError('Invalid PR number')
    if record.get('run_id') != str(run_id) or record.get('attempt') != int(attempt):
        raise ValueError('Artifact does not match the triggering workflow run')
    base = _text(record.get('base_sha'), 40)
    candidate = _text(record.get('candidate_sha'), 40)
    head = _text(record.get('head_sha'), 40)
    if not all(SHA.fullmatch(x) for x in (base, candidate, head)):
        raise ValueError('Invalid revision SHA')
    rows = record.get('rows')
    if not isinstance(rows, list) or not 0 < len(rows) <= MAX_ROWS:
        raise ValueError('Missing or excessive benchmark rows')
    names=set(); cleaned=[]
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError('Invalid benchmark row')
        name = _text(row.get('name'), 300)
        if name in names:
            raise ValueError('Duplicate benchmark')
        names.add(name)
        cleaned.append({
            'name': name,
            'base': _number(row.get('base')),
            'candidate': _number(row.get('candidate')),
            'delta': _number(row.get('delta')),
            'base_stddev': _number(row.get('base_stddev')),
            'candidate_stddev': _number(row.get('candidate_stddev')),
            'timing_status': _text(row.get('timing_status', ''), 80),
            'pair': _physics(row.get('pair')),
            'scheduled': _physics(row.get('scheduled')),
            'base_workload': row.get('base_workload'),
            'candidate_workload': row.get('candidate_workload'),
            'base_observables': row.get('base_observables'),
            'candidate_observables': row.get('candidate_observables'),
        })
        for key in ('base_workload', 'candidate_workload', 'base_observables', 'candidate_observables'):
            value=cleaned[-1][key]
            if value is not None and (not isinstance(value,dict) or len(value)>100):
                raise ValueError('Invalid workload or observable metadata')
            if value is not None:
                # Canonical JSON enforces finite values and simple serialisable structures.
                json.dumps(value, allow_nan=False)
    for field in ('base_versions', 'candidate_versions', 'base_compiler', 'candidate_compiler'):
        value=record.get(field) or {}
        if not isinstance(value, dict) or len(value)>60:
            raise ValueError('Invalid version or compiler metadata')
        json.dumps(value, allow_nan=False)
    return {
        'schema': 1, 'pr': pr, 'branch': _text(record.get('branch',''), 200),
        'run_id': str(run_id), 'attempt': int(attempt),
        'datetime': _text(record.get('datetime',''), 64),
        'base_sha': base, 'candidate_sha': candidate, 'head_sha': head,
        'machine': _text(record.get('machine',''), 200),
        'scheduled_reference': _text(record.get('scheduled_reference',''), 64),
        'base_versions': record.get('base_versions') or {},
        'candidate_versions': record.get('candidate_versions') or {},
        'base_compiler': record.get('base_compiler') or {},
        'candidate_compiler': record.get('candidate_compiler') or {},
        'rows': cleaned,
    }


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('artifact',type=Path)
    p.add_argument('--run-id',required=True,type=int)
    p.add_argument('--attempt',required=True,type=int)
    p.add_argument('--pr',required=True,type=int,
                   help='PR number verified independently from GitHub API')
    p.add_argument('--output',required=True,type=Path,
                   help='Trusted destination directory (not from the artifact)')
    args=p.parse_args(argv)
    if args.artifact.stat().st_size>MAX_BYTES:
        p.error('PR comparison artifact is too large')
    clean=validate(json.loads(args.artifact.read_text(encoding='utf-8')),
                   run_id=args.run_id,attempt=args.attempt)
    if clean['pr']!=args.pr:
        p.error('PR number does not match GitHub PR metadata')
    target=args.output/str(args.pr)/f'{args.run_id}-{args.attempt}.json'
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists():
        p.error('Refusing to overwrite an immutable PR comparison')
    target.write_text(json.dumps(clean,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    print(target)


if __name__=='__main__':
    main()