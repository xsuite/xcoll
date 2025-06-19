# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from pathlib import Path
import numpy as np
import pytest

import xtrack as xt
import xpart as xp
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts


num_part = 50000
num_turns = 3
path = Path(__file__).parent / 'data'


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("beam, plane", [
                            [1, 'H'],
                            [2, 'V'],
                            [1, 'V'],
                            [2, 'H']], ids=["B1H", "B2V", "B1V", "B2H"])
def test_impacts_from_line(beam, plane, test_context):
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    colldb = xc.CollimatorDatabase.from_yaml(path / 'colldb_lhc_run3.yaml', beam=beam)
    colldb.install_everest_collimators(verbose=True, line=line)
    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)

    impacts = xc.InteractionRecord.start(line=line, record_impacts=True, record_exits=True)
    line.build_tracker(_context=test_context)

    line.collimators.assign_optics()
    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    tw = line.twiss()
    part = line[tcp].generate_pencil(num_part, twiss=tw)

    line.scattering.enable()
    line.track(part, num_turns=num_turns, time=True, with_progress=1)
    line.scattering.disable()

    _assert_impacts(impacts, lengths=line.collimators.length)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_impacts_single_collimator(test_context):
    coll = xc.EverestCollimator(length=0.6, jaw=0.0013, material=xc.materials.MolybdenumGraphite,
                                emittance=3.5e-6, _context=test_context)

    x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
    px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

    impacts = xc.InteractionRecord.start(elements=[coll], names='TCP', record_impacts=True, record_exits=True)
    coll.track(part)
    part.sort(interleave_lost_particles=True)

    _assert_impacts(impacts, lengths=coll.length)


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
@pytest.mark.parametrize("R, side", [
                            [1, '+'],
                            [-1, '+'],
                            [1, '-'],
                            [-1, '-']], ids=["R>0 side=+ ", "R<0 side=+ ", "R>0 side=- ", "R<0 side=- "])
def test_impacts_single_crystal(R, side, test_context):
    coll = xc.EverestCrystal(length=0.002, material=xc.materials.SiliconCrystal, bending_angle=R*149e-6,
                        width=0.002, height=0.05, side=side, lattice='strip', jaw=0.001, _context=test_context)

    x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
    px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

    impacts = xc.InteractionRecord.start(elements=[coll], names='TCPCH', record_impacts=True, record_exits=True)
    coll.track(part)
    part.sort(interleave_lost_particles=True)

    _assert_impacts(impacts, expected_types=['Enter Jaw L', 'Exit Jaw'])


def _assert_impacts(impacts, expected_types=['Enter Jaw L', 'Enter Jaw R', 'Exit Jaw'], lengths=None):
    df = impacts.to_pandas()
    types = np.unique(df.interaction_type)
    assert np.all([type in expected_types for type in types])

    for this_type in ['Enter Jaw L', 'Enter Jaw R']:
        if this_type in expected_types:
            mask = df.interaction_type == this_type
            assert np.all(np.isclose(df.s_before[mask], 0.0, atol=1e-12) |
                        np.isclose(df.x_before[mask], 0.0, atol=1e-12))
            mask_all = mask & np.isclose(df.s_before, 0.0, atol=1e-12) & \
                       np.isclose(df.x_before, 0.0, atol=1e-12)
            assert mask_all.sum() < 3 # Allow maximally two particles at the corner
        else:
            assert this_type not in types

    if lengths:
        mask = df.interaction_type == 'Exit Jaw'
        if not isinstance(lengths, dict):
            lengths = {coll: lengths for coll in np.unique(df.collimator[mask])}
        assert np.all(np.isclose(df.s_before[mask], [lengths[coll] for coll in df.collimator[mask]], atol=1e-12) |
                    np.isclose(df.x_before[mask], 0.0, atol=1e-12))