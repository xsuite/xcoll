r"""Compare CPU physics of two Xcoll revisions using identical seeded inputs.

HOW TO USE
==========
The default matrix has 40 cases: eight standalone element workloads at
26, 450, and 7000 GeV/c, plus B1H/B1V/B2H/B2V with and without crystals at
450 and 7000 GeV/c. The LHC scan uses the same normalised Run 3 lattice at
both momenta; it is not an injection/flat-top optics comparison.

Generate a reference from a branch, tag, or commit (resolved to a full SHA):

    python benchmarks/everest_compare_physics.py reference \
        --ref main --output /tmp/physics-main

Use --ref v0.13.1 or --ref <full-sha> for an immutable revision. The command
creates and removes a temporary detached worktree. The output directory must
not exist. Both collections must use the same workload definitions, data, and
initial particles. Python, dependency and platform differences are reported,
not rejected. For a smaller reference use e.g.
--particles 20000 --lhc-particles 3000 --turns 5. To select a subset use
--cases crystal-channeling lhc --element-energies 450 7000
--lhc-energies 450 --beam 2 --plane H --crystals (only crystals).
--no-crystals selects only the non-crystal LHC cases; omission means both.

Manual worktrees (run from the candidate repository root):

    ROOT="$PWD"
    git worktree add --detach /tmp/xcoll-main main
    git worktree add --detach /tmp/xcoll-gpu EverestGPU
    SCRIPT="$ROOT/benchmarks/everest_compare_physics.py"
    DATA="$ROOT/tests/data"
    PYTHONPATH=/tmp/xcoll-main python "$SCRIPT" collect \
        --checkout /tmp/xcoll-main --data "$DATA" --output /tmp/physics-main
    PYTHONPATH=/tmp/xcoll-gpu python "$SCRIPT" collect \
        --checkout /tmp/xcoll-gpu --data "$DATA" --output /tmp/physics-gpu \
        --input-reference /tmp/physics-main
    python "$SCRIPT" compare /tmp/physics-main /tmp/physics-gpu
    git worktree remove /tmp/xcoll-main
    git worktree remove /tmp/xcoll-gpu

Automated test (current editable checkout, including uncommitted changes):

    XCOLL_EVEREST_PHYSICS_REFERENCE=/tmp/physics-main \
        python -m pytest tests/test_everest_physics_reference.py -v

Without an override the test reads tests/data/everest_physics_reference.
Use XCOLL_EVEREST_PHYSICS_REF_SHA to enforce an expected reference commit.
A manifest records the exact cases, Python, dependency versions and platform.
Candidates replay the stored initial particle coordinates and RNG states. The
comparison rejects incompatible workload parameters or initial inputs, but not
environment differences. --input-reference selects stored inputs for collect
and reference; particle counts, turns, interpolation and record capacity then
come from the selected reference cases. The PR workflow reads the accepted
reference and validation tools from a pinned main checkout, independently
of the candidate's tested merge/base SHA. CI never generates a reference.
XCOLL_EVEREST_PHYSICS_CHECKOUT selects the candidate source/data checkout
when the test harness lives in a different checkout.

For a permanent reference, generate from a validated commit (macOS or Linux),
inspect and review the results, then copy the entire reference directory
(including inputs/arrays.npz) into tests/data. Commit it
explicitly; never update an accepted reference merely to make a test pass.
Use a temporary directory or ignored location for ordinary generated output.
Compact existing schema-4 or schema-5 snapshots without rerunning physics:

    python benchmarks/everest_compare_physics.py compact /tmp/physics-v4 \
        --output /tmp/physics-v5

The comparison thresholds are screening criteria, not validated acceptance
criteria. Keep a higher-statistics reference for detailed studies.
"""
import argparse
from importlib import metadata
import json
from pathlib import Path
import subprocess
import sys
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
from everest_initial_particles import (
    initial_snapshot, initial_digest, restore_initial, read_initial,
)
from everest_snapshot_storage import (
    BANK_PATH, INPUT_STORAGE_VERSION, FINAL_FIELDS, SURVIVOR_FIELDS, compact_final, survivor_data,
    store_initial, load_initial,
)
from everest_interaction_statistics import (
    RECORD_COLUMNS,
    interaction_snapshot,
    compare_interactions,
)
from xcoll.interaction_record.interaction_record import (
    interaction_names,
)

WORKLOAD_VERSION = 5
DEFAULT_REFERENCE = DEFAULT_DATA / 'everest_physics_reference'
ELEMENT_ENERGIES = (26, 450, 7000)
LHC_ENERGIES = (450, 7000)
DEPENDENCIES = ('xobjects', 'xtrack', 'xpart', 'xdeps', 'numpy', 'scipy', 'pandas')
PROVENANCE_FIELDS = ('python', 'deps', 'platform')


def environment_metadata():
    deps = {p: metadata.version(p) for p in DEPENDENCIES}
    for package in ('xfields', 'xsuite'):
        try:
            deps[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            pass
    return dict(python=platform.python_version(), deps=deps,
                platform=dict(system=platform.system(), machine=platform.machine(),
                              release=platform.release()))


def environment_differences(ref, candidate):
    """Environment provenance is informational; physical inputs are not."""
    return [f'{key}: {ref.get(key)!r} -> {candidate.get(key)!r}'
            for key in PROVENANCE_FIELDS if ref.get(key) != candidate.get(key)]


def assert_compatible(ref, candidate):
    ignored = {'sha', *PROVENANCE_FIELDS}
    a = {k:v for k,v in ref.items() if k not in ignored}
    b = {k:v for k,v in candidate.items() if k not in ignored}
    if a != b:
        raise ValueError(f'Incompatible physical inputs for {ref.get("case")}:\n{a}\n{b}')


def lhc_start_names(line, initial):
    """Guard against silently tracking stored numeric indices in a changed line."""
    names = line.element_names
    indices = np.unique(initial['initial_at_element'])
    if len(indices) != 1 or not 0 <= indices[0] < len(names):
        raise ValueError('Expected a single valid LHC initial element index')
    start = int(initial['initial_start_tracking_at_element'])
    if not -1 <= start < len(names):
        raise ValueError('Invalid LHC tracking start index')
    return dict(at_element=names[int(indices[0])],
                start_tracking_at_element=names[start] if start >= 0 else None)



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


def load_manifest(path, *, version=WORKLOAD_VERSION):
    """Validate inventory and provenance before a test or comparison starts."""
    path = Path(path)
    manifest = json.loads((path / 'reference.json').read_text())
    cases = manifest.get('cases')
    if (manifest.get('workload_version') != version
            or not all(k in manifest for k in PROVENANCE_FIELDS)
            or not isinstance(cases, list) or not cases
            or any(not isinstance(c, str) or Path(c).name != c for c in cases)
            or len(set(cases)) != len(cases)
            or not re.fullmatch(r'[0-9a-f]{40}', str(manifest.get('sha', '')))):
        raise ValueError(f'Invalid reference manifest in {path}')
    actual = {p.stem for p in path.glob('*.npz')}
    if set(cases) != actual:
        raise ValueError(f'Incomplete reference {path}: {set(cases) ^ actual}')
    if version == WORKLOAD_VERSION and manifest.get("input_storage_version", 1) not in (1, INPUT_STORAGE_VERSION):
        raise ValueError(f'Unsupported input storage version in {path}')
    if version == WORKLOAD_VERSION and not (path / BANK_PATH).is_file():
        raise ValueError(f'Missing shared initial-particle bank in {path}')
    return manifest


def snapshot_info(path, case, manifest):
    with np.load(Path(path) / f'{case}.npz', allow_pickle=False) as z:
        info = json.loads(z['info'].item())
    for key in ('sha', *PROVENANCE_FIELDS, 'workload_version'):
        if info.get(key) != manifest.get(key):
            raise ValueError(f'{case}: {key} differs from manifest')
    family = info.get('family')
    if family not in (*ELEMENT_CASES, 'lhc'):
        raise ValueError(f'Unknown family in {case}')
    if not isinstance(info.get('n'), int) or info['n'] < 2:
        raise ValueError(f'Invalid particle count in {case}')
    if info.get('case') != case or case != case_id(
            family, info['energy_gev'], info.get('beam'), info.get('plane'),
            info.get('crystals')):
        raise ValueError(f'Inconsistent case metadata in {case}')
    return info


def snapshot_particles(p, n):
    ids = host(p, 'particle_id')
    keep = (ids >= 0) & (ids < n)
    order = np.argsort(ids[keep])
    rows = np.flatnonzero(keep)[order]
    if not np.array_equal(ids[rows], np.arange(n)):
        raise ValueError('Expected exactly one output row for every input primary')
    state = host(p, 'state')[rows]
    out = dict(state=state, at_element=host(p, 'at_element')[rows])
    survivor_rows = rows[state > 0]
    for field in SURVIVOR_FIELDS:
        out[f'survivor_{field}'] = host(p, field)[survivor_rows]
    return out


def save_case(out, case, p, n, initial, info, **extra):
    arrays = snapshot_particles(p, n)
    refs = store_initial(out, initial, n, rng_seed=info.get('particle_seed'))
    np.savez_compressed(out/f'{case}.npz', **arrays, **extra,
                        initial_refs=refs,
                        info=np.array(json.dumps(info, sort_keys=True)))


def particle_counts(args, specs):
    """Resolve default populations and explicit case overrides."""
    counts = {
        case: args.lhc_particles if spec['family'] == 'lhc' else args.particles
        for case, spec in specs}
    if any(n < 2 for n in counts.values()):
        raise ValueError('Physics particle populations must be >= 2')
    seen = set()
    for case, n in getattr(args, 'case_particles', []):
        if case not in counts:
            raise ValueError(f'Particle override refers to an unselected '
                             f'case: {case}')
        if case in seen:
            raise ValueError(f'Duplicate particle override: {case}')
        if n < 2:
            raise ValueError(f'Invalid particle count for {case}: {n}')
        counts[case] = n
        seen.add(case)
    return counts


def collect(args, *, allow_dirty=False, context=None,
            use_prebuilt=False):
    """Collect physics snapshots.

    By default, compile kernels from the imported Xcoll checkout.
    Prebuilt kernels may be used when they have explicitly been
    regenerated for that checkout.
    """
    with xo.settings.override(
        force_kernel_compilation=not use_prebuilt,
    ):
        return _collect(
            args,
            allow_dirty=allow_dirty,
            context=context,
        )

def _collect(args, *, allow_dirty=False, context=None):
    root = args.checkout.resolve()
    if not Path(xc.__file__).resolve().is_relative_to(root):
        raise SystemExit(f'Imported {xc.__file__}, expected Xcoll from {root}. '
                         'Set PYTHONPATH to the selected checkout.')
    def git(*a):
        return subprocess.check_output(['git','-C',str(root),*a], text=True).strip()
    if not allow_dirty:
        initial_status = git('status', '--porcelain')
        if initial_status:
            raise SystemExit('Physics comparison requires a clean checkout.')
    sha = git('rev-parse','HEAD')
    provenance = environment_metadata()
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    cpu = context if context is not None else xo.ContextCpu()
    specs = list(case_specs(args))
    counts = particle_counts(args, specs)
    input_reference = getattr(args, 'input_reference', None)
    if input_reference is not None:
        input_reference = Path(input_reference).resolve()
        input_manifest = load_manifest(input_reference)
        missing = set(counts) - set(input_manifest['cases'])
        if missing:
            raise ValueError(f'Input reference lacks selected cases: {sorted(missing)}')
    ids = [case for case, _ in specs]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate physics cases')
    if (out / 'reference.json').exists() or (out / BANK_PATH).exists() or any(out.glob('*.npz')):
        raise FileExistsError(f'Use a new output directory: {out}')
    for case, spec in specs:
        family = spec['family']
        energy = spec['energy_gev']
        print(f'Collecting {case} on {sha[:12]}', flush=True)
        source_info = None
        if input_reference is not None:
            source_info = snapshot_info(input_reference, case, input_manifest)
            with np.load(input_reference/f'{case}.npz', allow_pickle=False) as z:
                initial = load_initial(input_reference, z, source_info)
        n = source_info['n'] if source_info else counts[case]
        info = dict(case=case, n=n, sha=sha, **provenance,
                    workload_version=WORKLOAD_VERSION,
                    coordinate_seed=source_info['coordinate_seed'] if source_info else COORD_SEED,
                    particle_seed=source_info['particle_seed'] if source_info else PARTICLE_SEED,
                    **spec)
        if family == 'lhc':
            beam, plane, crystals = spec['beam'], spec['plane'], spec['crystals']
            line,tw,hashes = build_lhc(args.data, beam, crystals, cpu, energy_gev=energy)
            turns = source_info['turns'] if source_info else args.turns
            interpolation = source_info['interpolation'] if source_info else args.interpolation
            if source_info is None:
                p = make_lhc_pencil(line, tw, n, beam, plane)
                initial = initial_snapshot(p, n)
            p = restore_initial(initial, n, cpu)
            info.update(hashes=hashes, turns=turns, interpolation=interpolation,
                        start_names=lhc_start_names(line, initial),
                        input_sha256=initial_digest(initial, n),
                        p0c=float(initial['initial_p0c'][0]))
            if source_info:
                assert_compatible(source_info, info)
            line.xcoll.scattering.enable()
            reset_counters(line.elements)
            line.track(p,num_turns=turns)
            lm = make_lossmap(line,p,beam,interpolation)
            coll = lm.summary
            aper = lm.aperture_losses
            save_case(out,case,p,n,initial,info,
                coll_name=coll.name.to_numpy(dtype=str),coll_n=coll.n.to_numpy(),
                coll_e=coll.e.to_numpy(),aper_bin=np.asarray(aper['idx_bins']),
                aper_n=np.asarray(aper['n_bins']),aper_e=np.asarray(aper['e_bins']))
        else:
            el = make_element(family,cpu)
            el.name = 'element'  # Match the explicit InteractionRecord name.
            if source_info is None:
                p = make_particles(family,n,cpu,energy_gev=energy)
                initial = initial_snapshot(p, n)
            p = restore_initial(initial, n, cpu)
            info['input_sha256'] = initial_digest(initial, n)
            info['p0c'] = float(initial['initial_p0c'][0])
            # Record per-primary multiplicity and continuous longitudinal extent.
            info['events_per_particle'] = (source_info['events_per_particle']
                                           if source_info else args.events_per_particle)
            if source_info:
                assert_compatible(source_info, info)
            capacity = info['events_per_particle']*n
            # main still calls this argument capacity; EverestGPU calls it num_rows.
            size_key = ('num_rows' if 'num_rows' in inspect.signature(
                xc.InteractionRecord.__init__).parameters else 'capacity')
            rec = xc.InteractionRecord(
                elements=[el],
                names=['element'],
                **{size_key: capacity},
                columns=RECORD_COLUMNS,
                record_impacts=True,
                record_exits=True,
                record_scatterings=True,
            )
            try:
                el.track(p)
                nr = int(rec._index.num_recorded)
                if nr>=capacity:
                    raise ValueError('Record full; increase --events-per-particle')
                events = interaction_snapshot(rec, n=n, length=el.length, names=interaction_names)
                save_case(out, case, p, n, initial, info, **events,
                          deposited_energy=np.array(float(el._acc_ionisation_loss)))
            finally:
                rec.stop(elements=[el])
    if not allow_dirty:
        if (git('rev-parse', 'HEAD') != sha
                or git('status', '--porcelain') != initial_status):
            raise SystemExit('Checkout changed during collection; discard these snapshots.')
    manifest = dict(workload_version=WORKLOAD_VERSION, sha=sha,
                    **provenance, cases=ids)
    manifest['input_storage_version'] = INPUT_STORAGE_VERSION
    (out / 'reference.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    print(f'Saved {len(ids)} cases to {out}')


def mapping(z, k, v):
    return dict(zip(z[k].tolist(),z[v].tolist()))


def l1(a,b):
    return sum(abs(a.get(k,0)-b.get(k,0)) for k in a.keys()|b.keys())


def ks_screening_limit(n_ref, n_new, tolerance=0.02):
    # Approximate independent-sample allowance.
    return max(tolerance, 2.0 * np.sqrt(1/n_ref + 1/n_new))


def compare(args, *, cases=None):
    from scipy.stats import ks_2samp
    base, cand = args.base,args.candidate
    ma, mb = load_manifest(base), load_manifest(cand)
    for difference in environment_differences(ma, mb):
        print(f'Environment difference (informational): {difference}')
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
            assert_compatible(ia, ib)
            initial_a = load_initial(base, a, ia)
            initial_b = load_initial(cand, b, ib)
            n = ia['n']
            ida, sa = survivor_data(a, initial_a, n)
            idb, sb = survivor_data(b, initial_b, n)
            aa, bb = a['state'] > 0, b['state'] > 0
            _, common_a, common_b = np.intersect1d(
                ida, idb, assume_unique=True, return_indices=True)
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
            for field in (*SURVIVOR_FIELDS, 'kick_x'):
                x, y = sa[field], sb[field]
                if len(x) and len(y):
                    check(f'{field} survivor KS', ks_2samp(x, y).statistic, ks_screening_limit(len(x), len(y), args.ks_tol))
                    print(f'    mean {field}: {x.mean():.6g} / {y.mean():.6g}')
                    if len(common_a):
                        print(f'    same-ID max |difference|: {np.max(np.abs(x[common_a]-y[common_b])):.6g}')
                elif len(x)!=len(y):
                    check(f'{field} survivor KS',1.,args.ks_tol)
            if name.startswith('crystal-') and aa.any() and bb.any():
                check('channelled fraction (of survivors)',abs(
                    np.mean(sa['kick_x'] > 74.5e-6)-np.mean(sb['kick_x'] > 74.5e-6)),args.fraction_tol)
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
                compare_interactions(a, b, n, check, relative_tol=0.05, sigma=5.0)
                print(f'  deposited energy (eV): {a["deposited_energy"].item():.6g} / {b["deposited_energy"].item():.6g}')
                check('deposited energy difference / initial',abs(float(a['deposited_energy'])-float(b['deposited_energy']))/(n*ia['p0c']),args.fraction_tol)
    print(f'\n{failures} screening checks need review. Thresholds are not physics acceptance criteria.')
    return int(failures>0)


def compact_reference(source, output):
    """Migrate schema 4, or repack schema 5, without tracking."""
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f'Choose a new output directory: {output}')
    version = json.loads((source / 'reference.json').read_text())['workload_version']
    if version not in (4, WORKLOAD_VERSION):
        raise ValueError(f'Unsupported source workload version {version}')
    manifest = load_manifest(source, version=version)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='xcoll-compact-', dir=output.parent) as tmp:
        staging = Path(tmp) / 'snapshots'
        staging.mkdir()
        for case in manifest['cases']:
            info = snapshot_info(source, case, manifest)
            n = info['n']
            with np.load(source / f'{case}.npz', allow_pickle=False) as old:
                if version == 4:
                    initial = read_initial(old, info)
                    final = compact_final(old, n)
                else:
                    initial = load_initial(source, old, info)
                    final = {key: old[key] for key in ('state', 'at_element',
                             *(f'survivor_{f}' for f in SURVIVOR_FIELDS))}
                refs = store_initial(staging, initial, n, rng_seed=info.get('particle_seed'))
                new_info = dict(info, workload_version=WORKLOAD_VERSION)
                extra = {key: old[key] for key in old.files
                         if key not in FINAL_FIELDS and key != 'info'
                         and not key.startswith(('initial_', 'survivor_'))}
                dest = staging / f'{case}.npz'
                np.savez_compressed(dest, **final, **extra, initial_refs=refs,
                                    info=np.array(json.dumps(new_info, sort_keys=True)))
                # Verify every retained observable, including the derived kick.
                with np.load(dest, allow_pickle=False) as new:
                    recovered = load_initial(staging, new, new_info)
                    ids, values = survivor_data(new, recovered, n)
                    if version == 4:
                        original = {key: old[key][ids] for key in values}
                    else:
                        old_ids, original = survivor_data(old, initial, n)
                        if not np.array_equal(ids, old_ids):
                            raise ValueError(f'{case}: changed survivor IDs')
                    for key, value in values.items():
                        if not np.array_equal(value, original[key], equal_nan=True):
                            raise ValueError(f'{case}: changed survivor observable {key}')
                    for key in ('state', 'at_element', *extra):
                        if not np.array_equal(new[key], old[key],
                                              equal_nan=new[key].dtype.kind in 'fc'):
                            raise ValueError(f'{case}: changed observable {key}')
            print(f'Compacted {case}', flush=True)
        manifest['workload_version'] = WORKLOAD_VERSION
        manifest['input_storage_version'] = INPUT_STORAGE_VERSION
        (staging / 'reference.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
        load_manifest(staging)
        if output.exists():
            raise FileExistsError(f'Output destination now exists: {output}')
        staging.rename(output)
    print(f'Compacted {len(manifest["cases"])} cases to {output}')


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
            for case, count in args.case_particles:
                command += ['--case-particles', f'{case}={count}']
            if args.input_reference is not None:
                command += ['--input-reference', str(args.input_reference.resolve())]
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


def parse_case_particles(value):
    """Parse CASE=N, e.g. block-cfc-7000gev=500000."""
    try:
        case, count = value.rsplit('=', 1)
        count = int(count)
        if not case or count < 2:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError(
            'Expected CASE=INTEGER_AT_LEAST_2, e.g. '
            'block-cfc-7000gev=500000'
        ) from None
    return case, count


def add_collect_options(parser):
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--input-reference', type=Path,
                        help='Replay stored initial particles and per-case settings from this reference')
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
    parser.add_argument('--events-per-particle',type=int,default=12)
    parser.add_argument('--case-particles', action='append', type=parse_case_particles,
                        default=[], metavar='CASE=N',
                        help='Override the particle population for a specific '
                             'case. Can be specified multiple times.')


def main():
    p=argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
    sub=p.add_subparsers(dest='action',required=True)
    c=sub.add_parser('collect'); c.add_argument('--checkout',type=Path,required=True)
    add_collect_options(c)
    c.add_argument('--use-prebuilt', action='store_true', help='Use kernels already rebuilt for the selected checkout')
    r=sub.add_parser('reference', help='Generate reference from a Git ref in a temporary worktree')
    r.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    r.add_argument('--ref',required=True)
    add_collect_options(r)
    d=sub.add_parser('compare'); d.add_argument('base',type=Path);
    d.add_argument('candidate',type=Path)
    d.add_argument('--fraction-tol',type=float,default=.01);
    d.add_argument('--ks-tol',type=float,default=.02)
    d.add_argument('--map-tol',type=float,default=.05);
    m = sub.add_parser('compact', help='Migrate schema 4 or optimise schema 5 without retracking')
    m.add_argument('source', type=Path)
    m.add_argument('--output', type=Path, required=True)
    args=p.parse_args()
    if args.action == 'compact':
        compact_reference(args.source, args.output)
        return 0
    if args.action in ('collect', 'reference'):
        if min(args.particles,args.lhc_particles)<2 or min(args.turns,args.events_per_particle)<1 or not np.isfinite(args.interpolation) or args.interpolation<=0:
            p.error('particle counts must be >= 2; turns, capacity and interpolation must be positive')
        if args.action == 'collect':
            collect(args, use_prebuilt=args.use_prebuilt)
        else:
            reference(args)
        return 0
    if any(not np.isfinite(tol) or tol < 0 for tol in
           (args.fraction_tol,args.ks_tol,args.map_tol)):
        p.error('comparison tolerances must be finite and non-negative')
    return compare(args)


if __name__=='__main__':
    sys.exit(main())
