"""Serial CPU / OpenMP / CUDA crossover study with optional headless PNG.

Run from the repository root: python benchmarks/profile_gpu.py --help
All three contexts use identical seeded particles. Timings exclude initial
compilation and preparation, include synchronous tracking, and report both
tracking-only and copy + tracking + (for LHC) complete loss-map crossovers.
Tables are repeated at the end. If Matplotlib is installed, a PNG is saved
under benchmarks/results/ (or --plot PATH); --no-plot disables plotting.
No GUI or display server is required.
"""
import argparse
import contextlib
from datetime import datetime, timezone
import sys
from statistics import median
from time import perf_counter
import os
import platform
from pathlib import Path

import numpy as np
import xobjects as xo
import xcoll as xc
from everest_workloads import (DEFAULT_DATA, build_lhc, make_lhc_pencil,
    make_lossmap, make_element, make_particles, reset_counters, host)


def quiet_call(fn):
    # Xcoll/Xtrack can print progress even in a non-verbose loss-map call.
    with open(os.devnull, 'w') as sink, contextlib.redirect_stdout(sink):
        return fn()


def timed(fn, sync):
    sync()
    start = perf_counter()
    result = fn()
    sync()
    return perf_counter()-start, result


def omp_team_size(context):
    """Measure a real OpenMP team, not os.cpu_count() or a requested limit."""
    context.add_kernels(sources=[r'''
        #include <stdint.h>
        #include <omp.h>
        void xcoll_profile_omp_team(int32_t* size) {
            #pragma omp parallel
            {
                #pragma omp single
                *size = omp_get_num_threads();
            }
        }
    '''], kernels={'xcoll_profile_omp_team': xo.Kernel(
        args=[xo.Arg(xo.Int32, pointer=True, name='size')])}, specialize=False)
    size = np.zeros(1, dtype=np.int32)
    context.kernels.xcoll_profile_omp_team(size=size)
    return int(size[0]), int(context.omp_get_max_threads())


def run_case(initial, line_or_element, context, sync, rounds, turns=None,
             beam=1, interpolation=0.2):
    """Reset outside timing; measure input copy, track, and analysis separately."""
    is_lhc = turns is not None
    elements = ([el for el in line_or_element.elements
                 if hasattr(el, '_acc_ionisation_loss')] if is_lhc
                else [line_or_element])
    samples = []
    lost = None
    warmup = 0.
    for i in range(rounds+1):
        start = perf_counter()
        reset_counters(elements)
        t_copy, p = timed(lambda: initial.copy(_context=context), sync)
        track = (lambda: line_or_element.track(p, num_turns=turns)) if is_lhc else (
                 lambda: line_or_element.track(p))
        t_track, _ = timed(track, sync)
        t_map = 0.
        if is_lhc:
            t_map, lm = timed(lambda: make_lossmap(line_or_element, p, beam,
                                                   interpolation), sync)
            lost = lm.num_losses
        else:
            lost = np.count_nonzero(host(p, 'state') <= 0)
        if i:
            samples.append((t_copy, t_track, t_map, t_copy+t_track+t_map))
        else:
            warmup = perf_counter()-start
    return tuple(median(v) for v in zip(*samples)), int(lost), warmup


def crossover(rows, cpu, phase):
    """First sampled count with a GPU win at every remaining sampled count."""
    for i, (n, _) in enumerate(rows):
        if all(r[cpu][0][phase] > r['GPU'][0][phase] for _, r in rows[i:]):
            return n
    return None


CONTEXTS = ('CPU', 'OMP', 'GPU')


def table_phases(workload):
    phases = [('copy', 0), ('track', 1)]
    if workload == 'lhc':
        phases.append(('map', 2))
    return phases + [('total', 3)]


def print_table_header(workload):
    print(f'{"N":>9s}' + ''.join(f'{name+" "+phase:>13s}'
          for phase, _ in table_phases(workload) for name in CONTEXTS)
          + f'{"CPU/G track":>13s}{"OMP/G track":>13s}'
            f'{"CPU/G total":>13s}{"OMP/G total":>13s}', flush=True)


def print_table_row(workload, n, results):
    print(f'{n:9d}' + ''.join(f'{results[name][0][idx]:13.6f}'
          for _, idx in table_phases(workload) for name in CONTEXTS)
          + ''.join(f'{results[name][0][idx]/results["GPU"][0][idx]:13.2f}'
                    for idx in (1, 3) for name in ('CPU', 'OMP')), flush=True)
    print(' ' * 9 + 'losses CPU/OMP/GPU: ' + '/'.join(str(results[k][1]) for k in CONTEXTS)
          + '; excluded warm-up (s): ' + '/'.join(f'{results[k][2]:.3f}' for k in CONTEXTS),
          flush=True)


def print_crossovers(rows):
    for phase, label in [(1, 'tracking'), (3, 'input copy + tracking + analysis')]:
        for name in ('CPU', 'OMP'):
            n = crossover(rows, name, phase)
            print(f'{label} crossover GPU vs {name} (first sustained sampled win): '
                  + (str(n) if n is not None else 'not reached'))


def print_summary(reports, args, description):
    print('\n' + '=' * 78 + '\nFINAL GPU PROFILING SUMMARY\n' + '=' * 78)
    print(description)
    for workload, rows in reports.items():
        print(f'\n{workload.upper()} — median of {args.rounds} warmed runs')
        print_table_header(workload)
        for n, results in rows:
            print_table_row(workload, n, results)
        print_crossovers(rows)
    print('Times in seconds. Speed-ups > 1 favour CUDA. One-time preparation, warm-up, '
          'and pencil generation excluded. Total includes input copy, tracking, and '
          'loss-map analysis (including internal transfers and CPU refinement).')


def save_plot(reports, path, args, gpu_name, omp_threads):
    """Save the measured data without importing a GUI backend or showing a window."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print(f'PNG not generated: Matplotlib unavailable ({exc}). '
              'Install matplotlib to enable plotting.', file=sys.stderr)
        return False

    # One row per workload, with tracking, complete warmed job, and speed-up.
    fig, axes = plt.subplots(len(reports), 3, squeeze=False,
                             figsize=(15, 4.5*len(reports)), layout='constrained')
    labels = {'CPU': 'Serial CPU', 'OMP': f'OpenMP ({omp_threads} threads)',
              'GPU': 'CUDA'}
    styles = {'CPU': ('C0', 'o'), 'OMP': ('C1', 's'), 'GPU': ('C2', '^')}
    try:
        for i, (workload, rows) in enumerate(reports.items()):
            counts = [n for n, _ in rows]
            title = ('EverestBlock' if workload == 'block' else
                     f'LHC B{args.beam}{args.plane}, {args.turns} turns, '
                     f'{"with" if args.crystals else "without"} crystals')
            for j, (phase, heading) in enumerate(((1, 'Tracking'), (3, 'Total warmed job'))):
                ax = axes[i, j]
                for name in CONTEXTS:
                    colour, marker = styles[name]
                    ax.plot(counts, [r[name][0][phase] for _, r in rows],
                            color=colour, marker=marker, label=labels[name])
                ax.set_title(f'{title}\n{heading}')
                ax.set_ylabel('Median runtime [s]')
                ax.set_yscale('log')
                ax.legend(fontsize='small')
            ax = axes[i, 2]
            for name in ('CPU', 'OMP'):
                colour, marker = styles[name]
                for phase, suffix, linestyle in ((1, 'track', '--'), (3, 'total', '-')):
                    ax.plot(counts, [r[name][0][phase]/r['GPU'][0][phase]
                                     for _, r in rows], color=colour, marker=marker,
                            linestyle=linestyle, label=f'{labels[name]}/GPU {suffix}')
            ax.axhline(1., color='0.5', linewidth=1, linestyle=':')
            ax.set_title(f'{title}\nGPU speed-up')
            ax.set_ylabel('CPU time / GPU time')
            ax.set_yscale('log')
            ax.legend(fontsize='small')
            for ax in axes[i]:
                ax.set_xscale('log')
                ax.set_xlabel('Initial particles')
                ax.grid(True, which='both', alpha=0.25)
        fig.suptitle(f'Xcoll {xc.__version__} | {gpu_name} | '
                     f'{args.rounds} warmed rounds (median)')
        path = Path(path).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, format='png')
        print(f'Saved GPU profiling plot to {path}')
        return True
    finally:
        plt.close(fig)


def report(workload, counts, cpu, omp, gpu, cp, args):
    is_lhc = workload == 'lhc'
    print(f'\n{workload.upper()} — median of {args.rounds} warmed runs')
    sync_cpu = lambda: None
    sync_gpu = cp.cuda.runtime.deviceSynchronize
    contexts = {'CPU': cpu, 'OMP': omp, 'GPU': gpu}
    syncs = {'CPU': sync_cpu, 'OMP': sync_cpu, 'GPU': sync_gpu}
    start = perf_counter()
    if is_lhc:
        line, tw, _ = quiet_call(lambda: build_lhc(args.data, args.beam,
                                                   args.crystals, cpu))
        elements = {'CPU': line}
        for name in ('OMP', 'GPU'):
            ctx = contexts[name]
            other = quiet_call(lambda: line.copy(_context=ctx))
            quiet_call(lambda: other.build_tracker(_context=ctx))
            quiet_call(lambda: other.xcoll.scattering.enable())
            elements[name] = other
    else:
        elements = {name: make_element('block-cfc', ctx)
                    for name, ctx in contexts.items()}
    sync_gpu()
    print(f'One-time line/element preparation (all contexts): {perf_counter()-start:.2f} s '
          '(excluded from table)', flush=True)
    print_table_header(workload)
    rows = []
    for n in counts:
        initial = quiet_call(lambda: make_lhc_pencil(line, tw, n, args.beam,
                         args.plane)) if is_lhc else make_particles('block-cfc', n, cpu)
        # Identical initial bunch; rotate the context order between counts.
        results = {}
        order = list(contexts)
        offset = len(rows) % len(order)
        order = order[offset:] + order[:offset]
        for name in order:
            results[name] = quiet_call(lambda: run_case(initial, elements[name],
                 contexts[name], syncs[name], args.rounds,
                 args.turns if is_lhc else None, args.beam, args.interpolation))
        rows.append((n, results))
        print_table_row(workload, n, results)
    print_crossovers(rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workload', choices=['block','lhc','both'], default='both')
    parser.add_argument('--particles', type=int, nargs='+')
    parser.add_argument('--rounds', type=int, default=3)
    parser.add_argument('--turns', type=int, default=3)
    parser.add_argument('--cpu-threads', default='auto', metavar='auto|N',
                        help='OpenMP comparison thread count (default auto); serial is always included')
    parser.add_argument('--device', type=int, default=0)
    parser.add_argument('--beam', type=int, choices=[1,2], default=1)
    parser.add_argument('--plane', choices=['H','V'], default='H')
    parser.add_argument('--crystals', action='store_true')
    parser.add_argument('--interpolation', type=float, default=0.2)
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    plot_group = parser.add_mutually_exclusive_group()
    plot_group.add_argument('--plot', type=Path, metavar='PATH',
                            help='PNG output path (default benchmarks/results/gpu-profile-<UTC timestamp>.png)')
    plot_group.add_argument('--no-plot', action='store_true',
                            help='Do not generate a PNG')
    args = parser.parse_args()
    try:
        threads = 'auto' if args.cpu_threads == 'auto' else int(args.cpu_threads)
        if threads != 'auto' and threads < 1:
            raise ValueError
    except ValueError:
        parser.error('--cpu-threads must be auto or a positive integer')
    if min(args.rounds,args.turns)<1 or args.interpolation<=0 or (args.particles and min(args.particles)<1):
        parser.error('rounds, turns, interpolation and particle counts must be positive')
    import cupy as cp
    cp.cuda.Device(args.device).use()
    gpu = xo.ContextCupy(device=args.device)
    cpu = xo.ContextCpu(omp_num_threads=0)
    omp = xo.ContextCpu(omp_num_threads=threads)
    team, maximum = quiet_call(lambda: omp_team_size(omp))
    available = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count()
    name = cp.cuda.runtime.getDeviceProperties(args.device)['name']
    if isinstance(name, bytes):
        name = name.decode()
    description = (f'Xcoll {xc.__version__}; {platform.platform()}; serial CPU=1 thread; '
          f'OpenMP requested={threads}, observed team={team}, runtime maximum={maximum}, '
          f'available logical CPUs={available}; CUDA device={args.device} {name}; '
          f'turns={args.turns}; beam={args.beam}; plane={args.plane}; '
          f'crystals={args.crystals}; interpolation={args.interpolation} m')
    print(description)
    print('OpenMP team size is measured in a parallel region, not a count of distinct '
          'physical cores. Runtime affinity and dynamic teams can affect tracking.')
    reports = {}
    for workload in (['block','lhc'] if args.workload=='both' else [args.workload]):
        counts = args.particles or ([100,1000,10000,100000,1000000] if workload=='block'
                                    else [100,1000,10000,50000,100000])
        reports[workload] = report(workload, sorted(set(counts)), cpu, omp, gpu, cp, args)
    if not args.no_plot:
        stamp = datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S-%f')
        path = args.plot or Path(__file__).resolve().parent / 'results' / f'gpu-profile-{stamp}.png'
        try:
            save_plot(reports, path, args, name, team)
        except (OSError, RuntimeError, ValueError) as exc:
            print(f'PNG not generated: {exc}', file=sys.stderr)
    # Print the tables last, even if the optional plot could not be saved.
    print_summary(reports, args, description)


if __name__ == '__main__':
    main()
