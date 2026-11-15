"""Reproducible workloads shared by GPU benchmarks and standalone validation.

No pytest dependency: these helpers can be used against either Xcoll checkout.
"""
from pathlib import Path
import hashlib
import numpy as np
import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc

COORD_SEED = 20261109
PARTICLE_SEED = 13579
CRYSTAL_CASES = ('channeling', 'volume-reflection', 'amorphous', 'miss')
ELEMENT_CASES = ('block-cfc', 'block-tungsten', 'collimator-cfc',
                 'collimator-tungsten', *('crystal-'+c for c in CRYSTAL_CASES))
DEFAULT_DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'


def make_element(case, context):
    if case.startswith('crystal-'):
        return xc.EverestCrystal(length=0.002, material=xc.materials.Silicon,
            bending_angle=149e-6, width=0.002, height=0.05, side='+',
            jaw=0.001, miscut=0., lattice='strip', _context=context)
    material = (xc.materials.Tungsten if case.endswith('tungsten')
                else xc.materials.CarbonFibreComposite)
    if case.startswith('block-'):
        return xc.EverestBlock(length=0.6, material=material, _context=context)
    if case.startswith('collimator-'):
        return xc.EverestCollimator(length=0.6, material=material,
                                   jaw=[0.001, -0.001], _context=context)
    raise ValueError(f'Unknown element case: {case}')


def make_particles(case, n, context, energy_gev=None):
    rng = np.random.default_rng(COORD_SEED)
    if case.startswith('crystal-'):
        kind = case.removeprefix('crystal-')
        x = rng.normal(1.5e-3, 75e-6, n)
        if kind == 'channeling':
            px = rng.normal(0, 5e-6, n)
        elif kind == 'volume-reflection':
            px = rng.normal(75e-6, 15e-6, n)
        elif kind == 'amorphous':
            px = np.concatenate([rng.normal(180e-6, 5e-6, n//2),
                                 rng.normal(40e-6, 5e-6, n-n//2)])
        elif kind == 'miss':
            x = rng.normal(0, 0.1e-3, n)
            px = rng.normal(0, 30e-6, n)
        else:
            raise ValueError(case)
        p0c, y_sigma, py_sigma = 400e9, 1e-3, 5e-6
    else:
        x = (np.zeros(n) if case.startswith('block-')
             else rng.normal(1.025e-3, 8e-6, n))
        px = rng.normal(0, 3e-6, n)
        p0c, y_sigma, py_sigma = 450e9, 1e-4, 3e-6
    if energy_gev is not None:
        p0c = energy_gev * 1e9
    p = xp.Particles(p0c=p0c, x=x, px=px, y=rng.normal(0, y_sigma, n),
                     py=rng.normal(0, py_sigma, n), _context=context)
    seed_particles(p, n)
    return p


def seed_particles(p, n):
    p._init_random_number_generator(
        seeds=PARTICLE_SEED + np.arange(n, dtype=np.uint32))


def reset_counters(elements):
    for el in elements:
        if hasattr(el, '_acc_ionisation_loss'):
            el._acc_ionisation_loss = 0.
            el._acc_ionisation_loss_sec = 0.


def build_lhc(data=DEFAULT_DATA, beam=1, crystals=False, context=None,
              energy_gev=None):
    """Prepare real Run 3 optics on CPU. Reuse this line across particle counts."""
    data = Path(data)
    sequence = data / f'sequence_lhc_run3_b{beam}.json'
    database = data / 'colldb_lhc_run3_ir7.yaml'
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (sequence, database)}
    env = xt.load(sequence)
    line = env[f'lhcb{beam}']
    if energy_gev is not None:
        # Fixed normalised Run 3 lattice, with a new proton reference momentum.
        # This is an energy scan, not a separately matched injection lattice.
        line.particle_ref = xp.Particles(p0c=energy_gev * 1e9,
                                         mass0=xp.PROTON_MASS_EV)
    db = xc.CollimatorDatabase.from_yaml(database, beam=beam,
                                         ignore_crystals=not crystals)
    db.install_everest_collimators(line=line)
    line.xcoll.scattering.disable()
    line.build_tracker(_context=context or xo.ContextCpu())
    tw = line.twiss()
    line.xcoll.collimators.assign_optics(twiss=tw)
    if crystals:
        line.xcoll.collimators.align_to_beam_divergence()
    return line, tw, hashes


def make_lhc_pencil(line, tw, n, beam=1, plane='H'):
    tcp = f"tcp.{'c' if plane == 'H' else 'd'}6{'l' if beam == 1 else 'r'}7.b{beam}"
    line.xcoll.scattering.disable()
    old = np.random.get_state()
    try:
        np.random.seed(COORD_SEED)
        p = line.xcoll.collimators.generate_pencil(tcp, n, twiss=tw)
    finally:
        np.random.set_state(old)
    seed_particles(p, n)
    line.xcoll.scattering.enable()
    return p


def make_lossmap(line, p, beam=1, interpolation=0.2):
    return xc.LossMap(line=line, part=p, line_is_reversed=(beam == 2),
                      interpolation=interpolation, verbose=False)


def host(p, field):
    return np.asarray(p._context.nparray_from_context_array(getattr(p, field)))
