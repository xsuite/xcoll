r"""Compare CPU physics of two Xcoll revisions using identical seeded inputs.

HOW TO USE
==========
Run from the repository root, with the same pinned Python/Xsuite/NumPy/SciPy
and Pandas environment for both revisions. Use the same version of this script
and everest_workloads.py, and the same LHC test data for every collection.

Generate a reference directly from a branch, tag, or commit (resolved to a
full SHA before checkout), using a temporary detached Git worktree::

    python benchmarks/compare_everest_physics.py reference \
        --ref main --output /tmp/physics-main --crystals

Use --ref v0.13.1 or --ref <full-commit-sha> to pin the reference. The command
removes its temporary worktree, and never switches the current checkout.
The output directory must not already exist. Save the reported SHA with the
reference; do not silently regenerate it when main moves.

Alternatively, to collect both revisions explicitly in separate worktrees::

    ROOT="$PWD"
    git worktree add --detach /tmp/xcoll-main main
    git worktree add --detach /tmp/xcoll-gpu EverestGPU
    SCRIPT="$ROOT/benchmarks/compare_everest_physics.py"
    DATA="$ROOT/tests/data"
    PYTHONPATH=/tmp/xcoll-main python "$SCRIPT" collect \
        --checkout /tmp/xcoll-main --data "$DATA" \
        --output /tmp/physics-main --crystals
    PYTHONPATH=/tmp/xcoll-gpu python "$SCRIPT" collect \
        --checkout /tmp/xcoll-gpu --data "$DATA" \
        --output /tmp/physics-gpu --crystals
    python "$SCRIPT" compare /tmp/physics-main /tmp/physics-gpu
    git worktree remove /tmp/xcoll-main
    git worktree remove /tmp/xcoll-gpu

The script verifies which Xcoll is imported; the worktrees must contain the
clean revisions you intend to test. The workload definitions live beside this
script, not in the worktrees being measured. Use the same --cases, --particles,
--lhc-particles, --turns, --beam, --plane, --crystals and --interpolation flags
for both collections. For example, --cases lhc --beam 2 --plane H --crystals
selects the B2 horizontal LHC case. Use distinct directories for B1 and B2.

Automated regression test, using the generated reference and the current
editable checkout (local uncommitted candidate changes are allowed)::

    XCOLL_EVEREST_PHYSICS_REFERENCE=/tmp/physics-main \
        python -m pytest tests/test_everest_physics_reference.py -v

Optionally set XCOLL_EVEREST_PHYSICS_REF_SHA to the expected full reference
SHA. The test runs one parametrised case per reference snapshot, using the
reference's recorded workload parameters and current checkout's test data.
A missing reference configuration skips the opt-in test; an invalid supplied
reference, mismatched input/dependency, or failed comparison fails the test.

Each NPZ contains per-primary states/coordinates (sorted by particle ID),
interaction counts or loss-map bins, and input/provenance metadata. The
comparison rejects incompatible snapshots and returns a nonzero status for
screening failures. The configurable tolerances are preliminary screening
criteria, not a validated physics acceptance specification. NPZ files are
local reference artefacts; keep them outside the repository or in ignored
locations. Never regenerate an accepted reference from a failing candidate.
"""
import argparse
from importlib import metadata
import json
from pathlib import Path
import subprocess
import sys
import hashlib
import inspect
import os
import tempfile
import platform
import numpy as np
import xobjects as xo
import xcoll as xc
from everest_workloads import (DEFAULT_DATA, ELEMENT_CASES, COORD_SEED,
    PARTICLE_SEED, make_element, make_particles, build_lhc, make_lhc_pencil,
    make_lossmap, reset_counters, host)

FIELDS = ('x','px','y','py','zeta','delta','energy','s','state','at_element')


def digest(p):
    h = hashlib.sha256()
    for field in ('particle_id','x','px','y','py','zeta','delta','p0c','state','at_element',
                  '_rng_s1','_rng_s2','_rng_s3','_rng_s4'):
        h.update(np.ascontiguousarray(host(p,field)).tobytes())
    return h.hexdigest()


def snapshot_particles(p, n, initial_px):
    ids = host(p,'particle_id')
    keep = (ids >= 0) & (ids < n)
    order = np.argsort(ids[keep])
    ids = ids[keep][order]
    if not np.array_equal(ids, np.arange(n)):
        raise ValueError('Expected exactly one output row for every input primary')
    out = {f: host(p,f)[keep][order] for f in FIELDS}
    out['particle_id'] = ids
    out['kick_x'] = out['px']-initial_px[ids]
    return out


def save_case(out, case, p, n, initial_px, info, **extra):
    arrays = snapshot_particles(p,n,initial_px)
    np.savez_compressed(out/f'{case}.npz', **arrays, **extra,
                        info=np.array(json.dumps(info, sort_keys=True)))


def collect(args, *, allow_dirty=False):
    root = args.checkout.resolve()
    if not Path(xc.__file__).resolve().is_relative_to(root):
        raise SystemExit(f'Imported {xc.__file__}, expected Xcoll from {root}. '
                         'Set PYTHONPATH to the selected checkout.')
    def git(*a):
        return subprocess.check_output(['git','-C',str(root),*a], text=True).strip()
    initial_status = git('status','--porcelain')
    if initial_status and not allow_dirty:
        raise SystemExit('Physics comparison requires a clean checkout.')
    sha = git('rev-parse','HEAD')
    deps = {p: metadata.version(p) for p in ('xobjects','xtrack','xpart','xdeps','numpy','scipy','pandas')}
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    cpu = xo.ContextCpu()
    cases = args.cases or [*ELEMENT_CASES,'lhc']
    if len(set(cases)) != len(cases):
        raise ValueError('Duplicate physics cases')
    for case in cases:
        print(f'Collecting {case} on {sha[:12]}',flush=True)
        n = args.lhc_particles if case=='lhc' else args.particles
        info = dict(case=case,n=n,sha=sha,deps=deps,python=platform.python_version(),
                    workload_version=2,
                    coordinate_seed=COORD_SEED,
                    particle_seed=PARTICLE_SEED)
        if (out/f'{case}.npz').exists():
            raise FileExistsError(f'Use a new output directory: {out/f"{case}.npz"}')
        if case=='lhc':
            line,tw,hashes = build_lhc(args.data,args.beam,args.crystals,cpu)
            p = make_lhc_pencil(line,tw,n,args.beam,args.plane)
            info.update(hashes=hashes,beam=args.beam,plane=args.plane,
                        turns=args.turns,crystals=args.crystals,interpolation=args.interpolation)
            info['input_sha256'] = digest(p)
            info['p0c'] = float(host(p,'p0c')[0])
            initial_px = host(p,'px').copy()
            reset_counters(line.elements)
            line.track(p,num_turns=args.turns)
            lm = make_lossmap(line,p,args.beam,args.interpolation)
            coll = lm.summary
            aper = lm.aperture_losses
            save_case(out,case,p,n,initial_px,info,
                coll_name=coll.name.to_numpy(dtype=str),coll_n=coll.n.to_numpy(),
                coll_e=coll.e.to_numpy(),aper_bin=np.asarray(aper['idx_bins']),
                aper_n=np.asarray(aper['n_bins']),aper_e=np.asarray(aper['e_bins']))
        else:
            el = make_element(case,cpu)
            p = make_particles(case,n,cpu)
            info['input_sha256'] = digest(p)
            info['p0c'] = float(host(p,'p0c')[0])
            initial_px = host(p,'px').copy()
            # Minimal record columns: count each physical event, not the full coordinates.
            info['events_per_particle'] = args.events_per_particle
            capacity = args.events_per_particle*n
            # main still calls this argument capacity; EverestGPU calls it num_rows.
            size_key = ('num_rows' if 'num_rows' in inspect.signature(
                xc.InteractionRecord.__init__).parameters else 'capacity')
            rec = xc.InteractionRecord(elements=[el], names=['element'],
                **{size_key: capacity}, columns=[], record_impacts=True,
                record_exits=True, record_scatterings=True)
            try:
                el.track(p)
                nr = int(rec._index.num_recorded)
                if nr>=capacity:
                    raise ValueError('Record full; increase --events-per-particle')
                events = rec.to_pandas()['interaction_type'].value_counts()
                save_case(out,case,p,n,initial_px,info,
                    event_name=events.index.to_numpy(dtype=str),
                    event_n=events.to_numpy(),
                    deposited_energy=np.array(float(el._acc_ionisation_loss)))
            finally:
                rec.stop(elements=[el])
    if git('rev-parse','HEAD')!=sha or git('status','--porcelain') != initial_status:
        raise SystemExit('Checkout changed during collection; discard these snapshots.')
    print(f'Saved {len(cases)} cases to {out}')


def mapping(z, k, v):
    return dict(zip(z[k].tolist(),z[v].tolist()))


def l1(a,b):
    return sum(abs(a.get(k,0)-b.get(k,0)) for k in a.keys()|b.keys())


def compare(args, *, cases=None):
    from scipy.stats import ks_2samp
    base, cand = args.base,args.candidate
    a_names={p.name for p in base.glob('*.npz')}
    b_names={p.name for p in cand.glob('*.npz')}
    if cases is not None:
        selected = {f'{case}.npz' for case in cases}
        if not selected <= a_names or not selected <= b_names:
            raise ValueError(f'Missing requested cases: {selected - (a_names & b_names)}')
        a_names = b_names = selected
    if not a_names or a_names!=b_names:
        raise SystemExit(f'Missing/inconsistent cases: {a_names ^ b_names}')
    failures=0
    for name in sorted(a_names):
        with np.load(base/name,allow_pickle=False) as a, np.load(cand/name,allow_pickle=False) as b:
            ia,ib=json.loads(a['info'].item()),json.loads(b['info'].item())
            if {k:v for k,v in ia.items() if k!='sha'}!={k:v for k,v in ib.items() if k!='sha'}:
                raise SystemExit(f'Incompatible inputs or dependencies for {name}:\n{ia}\n{ib}')
            n=ia['n']; aa=a['state']>0; bb=b['state']>0
            print(f'\n{name}: {ia["sha"][:12]} → {ib["sha"][:12]} (N={n})')
            def check(label, value, limit):
                nonlocal failures
                ok=np.isfinite(value) and value<=limit
                failures+=int(not ok)
                print(f'  {label:40s} {value:12.6g}  limit={limit:g}  {"OK" if ok else "REVIEW"}')
            check('survival fraction difference',abs(aa.mean()-bb.mean()),args.fraction_tol)
            # Per-ID equality is diagnostic: an RNG path change may decorrelate histories.
            print(f'  identical particle states: {np.mean(a["state"]==b["state"]):.6%}')
            for code in np.union1d(a['state'],b['state']):
                check(f'state {code} fraction difference',abs(
                    np.mean(a['state']==code)-np.mean(b['state']==code)),args.fraction_tol)
            for field in ('x','px','y','py','delta','energy','kick_x'):
                x,y=a[field][aa],b[field][bb]
                if len(x) and len(y):
                    check(f'{field} survivor KS',ks_2samp(x,y).statistic,args.ks_tol)
                    print(f'    mean {field}: {x.mean():.6g} / {y.mean():.6g}')
                    common=aa&bb
                    if common.any():
                        print(f'    same-ID max |difference|: {np.max(np.abs(a[field][common]-b[field][common])):.6g}')
                elif len(x)!=len(y):
                    check(f'{field} survivor KS',1.,args.ks_tol)
            if name.startswith('crystal-') and aa.any() and bb.any():
                check('channelled fraction (of survivors)',abs(
                    np.mean(a['kick_x'][aa]>74.5e-6)-np.mean(b['kick_x'][bb]>74.5e-6)),args.fraction_tol)
            if name=='lhc.npz':
                for k,v,label in [('coll_name','coll_n','collimator loss L1 / N'),
                                  ('aper_bin','aper_n','aperture-bin loss L1 / N')]:
                    x,y=mapping(a,k,v),mapping(b,k,v)
                    check(label,l1(x,y)/n,args.map_tol)
                    for loc in sorted(x.keys()|y.keys(),key=lambda k: abs(x.get(k,0)-y.get(k,0)),reverse=True)[:5]:
                        print(f'    {loc}: {x.get(loc,0)} / {y.get(loc,0)}')
                for k in ('coll_e','aper_e'):
                    x,y=a[k].sum(),b[k].sum()
                    print(f'  {k} total (eV): {x:.6g} / {y:.6g}')
                    check(f'{k} energy difference / initial',abs(x-y)/(n*ia['p0c']),args.fraction_tol)
            else:
                x,y=mapping(a,'event_name','event_n'),mapping(b,'event_name','event_n')
                for k in sorted(x.keys()|y.keys()):
                    # Counts per incident proton, not probabilities (multiple events possible).
                    check(f'{k} events / primary',abs(x.get(k,0)-y.get(k,0))/n,args.event_tol)
                print(f'  deposited energy (eV): {a["deposited_energy"].item():.6g} / {b["deposited_energy"].item():.6g}')
                check('deposited energy difference / initial',abs(float(a['deposited_energy'])-float(b['deposited_energy']))/(n*ia['p0c']),args.fraction_tol)
    print(f'\n{failures} screening checks need review. Thresholds are not physics acceptance criteria.')
    return int(failures>0)


def reference(args):
    """Collect into a new directory using a temporary, detached worktree."""
    repo = args.repo.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f'Choose a new reference directory: {output}')
    output.parent.mkdir(parents=True, exist_ok=True)

    def git(*cmd):
        return subprocess.check_output(['git', '-C', str(repo), *cmd], text=True).strip()

    sha = git('rev-parse', '--verify', '--end-of-options', args.ref + '^{commit}')
    print(f'Reference {args.ref} resolves to {sha}', flush=True)
    with tempfile.TemporaryDirectory(prefix='xcoll-physics-', dir=output.parent) as tmp:
        worktree = Path(tmp)/'checkout'
        staging = Path(tmp)/'snapshots'
        git('worktree', 'add', '--detach', str(worktree), sha)
        try:
            command = [sys.executable, str(Path(__file__).resolve()), 'collect',
                       '--checkout', str(worktree), '--output', str(staging),
                       '--data', str(args.data.resolve())]
            for option in ('particles', 'lhc_particles', 'turns', 'beam', 'plane',
                           'interpolation', 'events_per_particle'):
                command += ['--'+option.replace('_','-'), str(getattr(args, option))]
            if args.cases:
                command += ['--cases', *args.cases]
            if args.crystals:
                command += ['--crystals']
            env = os.environ.copy()
            env['PYTHONPATH'] = str(worktree) + (os.pathsep + env['PYTHONPATH']
                                                 if env.get('PYTHONPATH') else '')
            subprocess.run(command, env=env, cwd=worktree, check=True)
            if output.exists():
                raise FileExistsError(f'Reference destination now exists: {output}')
            staging.rename(output)
        finally:
            git('worktree', 'remove', '--force', str(worktree))
    print(f'Reference {sha} saved in {output}')


def add_collect_options(parser):
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--data',type=Path,default=DEFAULT_DATA)
    parser.add_argument('--cases',nargs='+',choices=[*ELEMENT_CASES,'lhc'])
    parser.add_argument('--particles',type=int,default=100_000)
    parser.add_argument('--lhc-particles',type=int,default=12_000)
    parser.add_argument('--turns',type=int,default=10)
    parser.add_argument('--beam',type=int,choices=[1,2],default=1)
    parser.add_argument('--plane',choices=['H','V'],default='H')
    parser.add_argument('--crystals',action='store_true')
    parser.add_argument('--interpolation',type=float,default=0.2)
    parser.add_argument('--events-per-particle',type=int,default=40)


def main():
    p=argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
    sub=p.add_subparsers(dest='action',required=True)
    c=sub.add_parser('collect'); c.add_argument('--checkout',type=Path,required=True)
    add_collect_options(c)
    r=sub.add_parser('reference', help='Generate reference from a Git ref in a temporary worktree')
    r.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    r.add_argument('--ref',required=True)
    add_collect_options(r)
    d=sub.add_parser('compare'); d.add_argument('base',type=Path); d.add_argument('candidate',type=Path)
    d.add_argument('--fraction-tol',type=float,default=.01); d.add_argument('--ks-tol',type=float,default=.02)
    d.add_argument('--map-tol',type=float,default=.05); d.add_argument('--event-tol',type=float,default=.03)
    args=p.parse_args()
    if args.action in ('collect', 'reference'):
        if min(args.particles,args.lhc_particles,args.turns,args.events_per_particle)<1 or args.interpolation<=0:
            p.error('particle counts, turns, capacity and interpolation must be positive')
        if args.action == 'collect':
            collect(args)
        else:
            reference(args)
        return 0
    if any(not np.isfinite(tol) or tol < 0 for tol in
           (args.fraction_tol,args.ks_tol,args.map_tol,args.event_tol)):
        p.error('comparison tolerances must be finite and non-negative')
    return compare(args)


if __name__=='__main__':
    sys.exit(main())
