"""Inspect recorded scattering statistics in Xcoll physics reference snapshots.

Usage from the Xcoll root:
    python inspect_everest_statistics.py tests/data/everest_physics_reference

Counts are recorded events/segments, not necessarily distinct particles.
The projected statistical resolution assumes independent Poisson counts; it is
only a planning approximation, particularly for continuous processes.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

PROCESSES = (
    'Absorbed', 'Multiple Coulomb Scattering', 'PN Elastic', 'PP Elastic',
    'Single Diffractive', 'Coulomb', 'Channelling', 'Dechannelling',
    'Volume Reflection Trans CH', 'Volume Reflection',
    'Volume Reflection Trans MCS', 'Multiple Coulomb Trans VR', 'Volume Capture',
)
CONTINUOUS = {'multiple_coulomb_scattering', 'channelling', 'multiple_coulomb_trans_vr'}


def key(name):
    return '_'.join(str(name).lower().replace('_', ' ').split())


def read_cases(path):
    cases = []
    for file in sorted(path.glob('*.npz')):
        with np.load(file, allow_pickle=False) as z:
            if 'event_name' not in z:
                continue  # LHC snapshots do not contain interaction counts.
            info = json.loads(z['info'].item())
            n = int(info['n'])
            if n < 1:
                raise ValueError(f'{file}: invalid particle count')
            counts = {}
            for name, count in zip(z['event_name'], z['event_n']):
                k = key(name)
                if k.startswith(('enter_jaw', 'exit_jaw')):
                    continue
                counts[k] = counts.get(k, 0) + int(count)
            cases.append((file.stem, n, counts))
    if not cases:
        raise ValueError(f'No element interaction snapshots found in {path}')
    return cases


def report(cases, target, sigma, proposed_n):
    labels = {key(p): p for p in PROCESSES}
    # Keep new/unknown interaction types rather than silently dropping them.
    for _, _, counts in cases:
        for k in counts:
            labels.setdefault(k, k.replace('_', ' '))
    print('Recorded counts by case (0 = not observed, not necessarily impossible)')
    print(f'{"Case":36s} {"N":>8s} ' + ' '.join(f'{i+1:>7d}' for i in range(len(labels))))
    for case, n, counts in cases:
        print(f'{case:36s} {n:8d} ' + ' '.join(f'{counts.get(k,0):7d}' for k in labels))
    print('\nColumn key:')
    for i, (k, name) in enumerate(labels.items(), 1):
        print(f'{i:2d} {name}' + (' (continuous segments)' if k in CONTINUOUS else ''))
    print('\nBest observed yield for each process:')
    print(f'{"Process":32s} {"Best case":36s} {"Events":>9s} {"Events/N":>10s} {"Projected":>12s} {"N target":>12s}')
    required_events = 2 * (sigma / target)**2
    for k, name in labels.items():
        available = [(case, n, counts.get(k,0)) for case,n,counts in cases if counts.get(k,0)>0]
        if not available:
            print(f'{name:32s} {"not observed in any case":36s}')
            continue
        case,n,c = max(available, key=lambda t:t[2]/t[1])
        rate = c/n
        required = math.ceil(required_events/rate)
        print(f'{name:32s} {case:36s} {c:9d} {rate:10.4g} {rate*proposed_n:12.0f} {required:12d}')
    print(f'\nProjected = estimated recorded events at N={proposed_n:,}.')
    print(f'N target = approximate incident count for {sigma:g}-sigma resolution of a {100*target:g}% relative change, assuming independent Poisson counts.')
    print('For continuous segments or repeated interactions per proton, use per-primary variances instead of Poisson statistics.')
    print('A rare process may require a deliberately targeted geometry or incidence distribution, rather than more particles in the current case.')


def main():
    p=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('reference', type=Path)
    p.add_argument('--relative-change', type=float, default=.10)
    p.add_argument('--sigma', type=float, default=5.)
    p.add_argument('--particles', type=int, default=500_000)
    a=p.parse_args()
    if not 0<a.relative_change<1 or not math.isfinite(a.sigma) or a.sigma<=0 or a.particles<1:
        p.error('relative-change must be between 0 and 1; sigma and particles must be positive')
    report(read_cases(a.reference),a.relative_change,a.sigma,a.particles)


if __name__=='__main__':
    main()
