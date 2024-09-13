# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pytest
from pathlib import Path

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts


num_part = 50000
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
    coll_manager = xc.CollimatorDatabase.from_yaml(path / f'colldb_lhc_run3.yaml', beam=beam)
    coll_manager.install_everest_collimators(verbose=True, line=line)
    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)

    impacts = xc.InteractionRecord.start(line=line, record_impacts=True, record_exits=True)
    line.build_tracker(_context=test_context)

    xc.assign_optics_to_collimators(line=line)
    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    tw = line.twiss()
    part = xc.generate_pencil_on_collimator(line, tcp, num_particles=num_part, tw=tw)

    xc.enable_scattering(line)
    line.track(part, num_turns=20, time=True, with_progress=1)
    xc.disable_scattering(line)

    df = impacts.to_pandas()
    types = np.unique(df.interaction_type)
    assert np.all([type in ['Enter Jaw L', 'Enter Jaw R', 'Exit Jaw'] for type in types])

    mask = df.interaction_type == 'Enter Jaw L'
    assert np.all(np.isclose(df.s_before[mask], 0.0, atol=1e-12) |
                  np.isclose(df.x_before[mask], 0.0, atol=1e-12))
    mask = df.interaction_type == 'Enter Jaw R'
    assert np.all(np.isclose(df.s_before[mask], 0.0, atol=1e-12) |
                  np.isclose(df.x_before[mask], 0.0, atol=1e-12))
    mask = df.interaction_type == 'Exit Jaw'
    assert np.all(np.isclose(df.s_before[mask], [line[coll].length for coll in df.collimator[mask]], atol=1e-12) |
                  np.isclose(df.x_before[mask], 0.0, atol=1e-12))


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

    impacts_coll = xc.InteractionRecord.start(elements=[coll], names='TCP', record_impacts=True, record_exits=True)
    coll.track(part)
    part.sort(interleave_lost_particles=True)

    df = impacts_coll.to_pandas()
    types = np.unique(df.interaction_type)
    assert np.all([type in ['Enter Jaw L', 'Enter Jaw R', 'Exit Jaw'] for type in types])

    mask = df.interaction_type == 'Enter Jaw L'
    assert np.all(np.isclose(df.s_before[mask], 0.0, atol=1e-12) |
                  np.isclose(df.x_before[mask], 0.0, atol=1e-12))
    mask = df.interaction_type == 'Enter Jaw R'
    assert np.all(np.isclose(df.s_before[mask], 0.0, atol=1e-12) |
                  np.isclose(df.x_before[mask], 0.0, atol=1e-12))
    mask = df.interaction_type == 'Exit Jaw'
    assert np.all(np.isclose(df.s_before[mask], coll.length, atol=1e-12) |
                  np.isclose(df.x_before[mask], 0.0, atol=1e-12))


@for_all_test_contexts(
    excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
)
def test_impacts_single_crystal(test_context):
    coll = xc.EverestCrystal(length=0.002, material=xc.materials.SiliconCrystal, bending_angle=149e-6,
                        width=0.002, height=0.05, side='+', lattice='strip', jaw=0.001, _context=test_context)

    x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
    px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

    impacts_cry = xc.InteractionRecord.start(elements=[coll], names='TCPCH', record_impacts=True, record_exits=True)

    coll.track(part)
    part.sort(interleave_lost_particles=True)

    df = impacts_cry.to_pandas()
    types = np.unique(df.interaction_type)
    assert np.all([type in ['Enter Jaw L', 'Exit Jaw'] for type in types])
    assert 'Enter Jaw R' not in types

    mask = df.interaction_type == 'Enter Jaw L'
    assert np.all(np.isclose(df.s_before[mask], 0.0, atol=1e-12) |
                  np.isclose(df.x_before[mask], 0.0, atol=1e-12))
    mask = df.interaction_type == 'Exit Jaw'
    assert np.all(np.isclose(df.s_before[mask], coll.length, atol=1e-12) |
                  np.isclose(df.x_before[mask], 0.0, atol=1e-12))
