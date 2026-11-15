r"""Compare CPU physics of two Xcoll revisions using identical seeded inputs.

HOW TO USE
==========
The default matrix has 40 cases: eight standalone element workloads at
26, 450, and 7000 GeV/c, plus B1H/B1V/B2H/B2V with and without crystals at
450 and 7000 GeV/c. The LHC scan uses the same normalised Run 3 lattice at
both momenta; it is not an injection/flat-top optics comparison.

Generate a reference from a branch, tag, or commit (resolved to a full SHA):

    python benchmarks/compare_everest_physics.py reference \
        --ref main --output /tmp/physics-main

Use --ref v0.13.1 or --ref <full-sha> for an immutable revision. The command
creates and removes a temporary detached worktree. The output directory must
not exist. Both collections must use the same workload definitions, data,
Python and dependency versions. For a smaller reference use e.g.
--particles 20000 --lhc-particles 3000 --turns 5. To select a subset use
--cases crystal-channeling lhc --element-energies 450 7000
--lhc-energies 450 --beam 2 --plane H --crystals (only crystals).
--no-crystals selects only the non-crystal LHC cases; omission means both.

Manual worktrees (run from the candidate repository root):

    ROOT="$PWD"
    git worktree add --detach /tmp/xcoll-main main
    git worktree add --detach /tmp/xcoll-gpu EverestGPU
    SCRIPT="$ROOT/benchmarks/compare_everest_physics.py"
    DATA="$ROOT/tests/data"
    PYTHONPATH=/tmp/xcoll-main python "$SCRIPT" collect \
        --checkout /tmp/xcoll-main --data "$DATA" --output /tmp/physics-main
    PYTHONPATH=/tmp/xcoll-gpu python "$SCRIPT" collect \
        --checkout /tmp/xcoll-gpu --data "$DATA" --output /tmp/physics-gpu
    python "$SCRIPT" compare /tmp/physics-main /tmp/physics-gpu
    git worktree remove /tmp/xcoll-main
    git worktree remove /tmp/xcoll-gpu

Automated test (current editable checkout, including uncommitted changes):

    XCOLL_EVEREST_PHYSICS_REFERENCE=/tmp/physics-main \
        python -m pytest tests/test_everest_physics_reference.py -v

Without an override the test reads tests/data/everest_physics_reference.
Use XCOLL_EVEREST_PHYSICS_REF_SHA to enforce an expected reference commit.
A manifest records the exact cases, Python and dependency versions. The test
rejects incomplete or incompatible snapshots. The PR workflow reads the
committed reference from the base checkout, never from the candidate; before
one exists it collects from the exact PR base SHA in a fresh worktree.

For a permanent reference, generate from a validated commit in the CI Python
and dependency environment, inspect and review the results, then copy the
reference directory (manifest and all NPZ files) into tests/data. Commit it
explicitly; never update an accepted reference merely to make a test pass.
Use a temporary directory or ignored location for ordinary generated output.
The comparison thresholds are screening criteria, not validated acceptance
criteria. Keep a higher-statistics reference for detailed studies.
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
import re
from itertools import product
import numpy as np
import xobjects as xo
import xcoll as xc
from everest_workloads import (DEFAULT_DATA, ELEMENT_CASES, COORD_SEED,
    PARTICLE_SEED, make_element, make_particles, build_lhc, make_lhc_pencil,
    make_lossmap, reset_counters, host)

FIELDS = ('x','px','y','py','zeta','delta','energy','s','state','at_element')
WORKLOAD_VERSION = 3
DEFAULT_REFERENCE = DEFAULT_DATA / 'everest_physics_reference'
ELEMENT_ENERGIES = (26, 450, 7000)
LHC_ENERGIES = (450, 7000)
DEPENDENCIES = ('xobjects', 'xtrack', 'xpart', 'xdeps', 'numpy', 'scipy', 'pandas')


def case_id(family, energy_gev, beam=None, plane=None, crystals=None):
    if family == 'lhc':
        mode = 'crystals' if crystals else 'no-crystals'
        return f'lhc-b{beam}-{plane.lower()}-{mode}-{energy_gev}gev'
    return f'{family}-{energy_gev}gev'


def case_specs(args):
    """Enumerate unique, deterministic workload IDs and their parameters."""
    families = args.cases or [*ELEMENT_CASES, 'lhc']
    if len(set(families)) != len(families):
        raise ValueError('Duplicate physics case families')
    for family in families:
        if family == 'lhc':
            for energy, beam, plane, crystals in product(
                    args.lhc_energies, args.beam, args.plane, args.crystal_modes):
                yield case_id(family, energy, beam, plane, crystals), dict(
                    family=family, energy_gev=energy, beam=beam, plane=plane,
                    crystals=crystals)
        else:
            for energy in args.element_energies:
                yield case_id(family, energy), dict(family=family, energy_gev=energy)


def load_manifest(path):
    """Validate inventory and provenance before a test or comparison starts."""
    path = Path(path)
    manifest = json.loads((path / 'reference.json').read_text())
    cases = manifest.get('cases')
    if (manifest.get('workload_version') != WORKLOAD_VERSION
            or not isinstance(cases, list) or not cases
            or any(not isinstance(c, str) or Path(c).name != c for c in cases)
            or len(set(cases)) != len(cases)
            or not re.fullmatch(r'[0-9a-f]{40}', str(manifest.get('sha', '')))):
        raise ValueError(f'Invalid reference manifest in {path}')
    actual = {p.stem for p in path.glob('*.npz')}
    if set(cases) != actual:
        raise ValueError(f'Incomplete reference {path}: {set(cases) ^ actual}')
    return manifest


def snapshot_info(path, case, manifest):
    with np.load(Path(path) / f'{case}.npz', allow_pickle=False) as z:
        info = json.loads(z['info'].item())
    for key in ('sha', 'python', 'deps', 'workload_version'):
        if info.get(key) != manifest.get(key):
            raise ValueError(f'{case}: {key} differs from manifest')
    family = info.get('family')
    if family not in (*ELEMENT_CASES, 'lhc'):
        raise ValueError(f'Unknown family in {case}')
    if info.get('case') != case or case != case_id(
            family, info['energy_gev'], info.get('beam'), info.get('plane'),
            info.get('crystals')):
        raise ValueError(f'Inconsistent case metadata in {case}')
    return info


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
    # Scope the compilation policy to this collection (pytest may run other tests).
    with xo.settings.override(force_kernel_compilation=True):
        return _collect(args, allow_dirty=allow_dirty)


def _collect(args, *, allow_dirty=False):
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
    deps = {p: metadata.version(p) for p in DEPENDENCIES}
    # Record optional packages when installed, too.
    for package in ('xfields', 'xsuite'):
        try:
            deps[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            pass
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    cpu = xo.ContextCpu()
    specs = list(case_specs(args))
    ids = [case for case, _ in specs]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate physics cases')
    if (out / 'reference.json').exists() or any(out.glob('*.npz')):
        raise FileExistsError(f'Use a new output directory: {out}')
    for case, spec in specs:
        family = spec['family']
        energy = spec['energy_gev']
        print(f'Collecting {case} on {sha[:12]}', flush=True)
        n = args.lhc_particles if family == 'lhc' else args.particles
        info = dict(case=case, n=n, sha=sha, deps=deps,
                    python=platform.python_version(), workload_version=WORKLOAD_VERSION,
                    coordinate_seed=COORD_SEED, particle_seed=PARTICLE_SEED, **spec)
        if family == 'lhc':
            beam, plane, crystals = spec['beam'], spec['plane'], spec['crystals']
            line,tw,hashes = build_lhc(args.data, beam, crystals, cpu, energy_gev=energy)
            p = make_lhc_pencil(line, tw, n, beam, plane)
            info.update(hashes=hashes, turns=args.turns, interpolation=args.interpolation)
            info['input_sha256'] = digest(p)
            info['p0c'] = float(host(p,'p0c')[0])
            initial_px = host(p,'px').copy()
            reset_counters(line.elements)
            line.track(p,num_turns=args.turns)
            lm = make_lossmap(line,p,beam,args.interpolation)
            coll = lm.summary
            aper = lm.aperture_losses
            save_case(out,case,p,n,initial_px,info,
                coll_name=coll.name.to_numpy(dtype=str),coll_n=coll.n.to_numpy(),
                coll_e=coll.e.to_numpy(),aper_bin=np.asarray(aper['idx_bins']),
                aper_n=np.asarray(aper['n_bins']),aper_e=np.asarray(aper['e_bins']))
        else:
            el = make_element(family,cpu)
            el.name = 'element'  # Match the explicit InteractionRecord name.
            p = make_particles(family,n,cpu,energy_gev=energy)
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
    manifest = dict(workload_version=WORKLOAD_VERSION, sha=sha,
                    python=platform.python_version(), deps=deps, cases=ids)
    (out / 'reference.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    print(f'Saved {len(ids)} cases to {out}')


def mapping(z, k, v):
    return dict(zip(z[k].tolist(),z[v].tolist()))


def l1(a,b):
    return sum(abs(a.get(k,0)-b.get(k,0)) for k in a.keys()|b.keys())


def compare(args, *, cases=None):
    from scipy.stats import ks_2samp
    base, cand = args.base,args.candidate
    ma, mb = load_manifest(base), load_manifest(cand)
    a_names={f'{c}.npz' for c in ma['cases']}
    b_names={f'{c}.npz' for c in mb['cases']}
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
            ia = snapshot_info(base, Path(name).stem, ma)
            ib = snapshot_info(cand, Path(name).stem, mb)
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
            if ia['family'] == 'lhc':
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
            for option in ('particles', 'lhc_particles', 'turns',
                           'interpolation', 'events_per_particle'):
                command += ['--'+option.replace('_','-'), str(getattr(args, option))]
            if args.cases:
                command += ['--cases', *args.cases]
            for option in ('element_energies', 'lhc_energies', 'beam', 'plane'):
                command += ['--'+option.replace('_','-'), *map(str, getattr(args, option))]
            if args.crystal_modes == [True]:
                command += ['--crystals']
            elif args.crystal_modes == [False]:
                command += ['--no-crystals']
            env = os.environ.copy()
            # A prebuilt kernel from the active environment might contain a
            # different Xcoll revision than this detached worktree.
            env['XSUITE_FORCE_KERNEL_COMPILATION'] = '1'
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
    parser.add_argument('--element-energies', nargs='+', type=int, choices=ELEMENT_ENERGIES, default=list(ELEMENT_ENERGIES))
    parser.add_argument('--lhc-energies', nargs='+', type=int, choices=LHC_ENERGIES, default=list(LHC_ENERGIES))
    parser.add_argument('--beam', nargs='+', type=int, choices=[1,2], default=[1,2])
    parser.add_argument('--plane', nargs='+', choices=['H','V'], default=['H','V'])
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--crystals', dest='crystal_modes', action='store_const', const=[True])
    modes.add_argument('--no-crystals', dest='crystal_modes', action='store_const', const=[False])
    parser.set_defaults(crystal_modes=[False, True])
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
        if min(args.particles,args.lhc_particles,args.turns,args.events_per_particle)<1 or not np.isfinite(args.interpolation) or args.interpolation<=0:
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
