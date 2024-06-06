# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest 
from pathlib import Path

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc
from xobjects.test_helpers import for_all_test_contexts
# ---------------------------------------------------------

# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
path = Path(__file__).parent / 'data'
@pytest.mark.parametrize("beam, plane", [
                            [1, 'H'],
                            [2, 'V'],
                            [1, 'V'],
                            [2, 'H']], ids=["B1H", "B2V", "B1V", "B2H"])

def test_impacts_coll(beam, plane):
    # Collimator
    line = xt.Line.from_json(path / f'sequence_lhc_run3_b{beam}.json')
    coll_manager = xc.CollimatorDatabase.from_yaml(path / f'colldb_lhc_run3.yaml', beam=beam)
    coll_manager.install_everest_collimators(verbose=True, line=line)
    df_with_coll = line.check_aperture()
    assert not np.any(df_with_coll.has_aperture_problem)

    impacts = xc.InteractionRecord.start(line)
    line.build_tracker()

    xc.assign_optics_to_collimators(line=line)
    tcp  = f"tcp.{'c' if plane=='H' else 'd'}6{'l' if beam==1 else 'r'}7.b{beam}"
    tw = line.twiss()
    part = xc.generate_pencil_on_collimator(line, tcp, 50000, tw=tw)

    line.discard_tracker()
    line.build_tracker(_context=xo.ContextCpu(omp_num_threads=12))

    xc.enable_scattering(line)
    line.track(part, num_turns=20, time=True, with_progress=1)
    xc.disable_scattering(line)

    line.discard_tracker()

    df = impacts.to_pandas()

    assert (np.allclose(df[df.interaction_type == 'Enter Jaw L']['ds'], 0.0, atol=1e-12) or 
            np.allclose(df[df.interaction_type == 'Enter Jaw L']['parent_x'], 0.0, atol=1e-12))
    assert (np.allclose(df[df.interaction_type == 'Enter Jaw R']['ds'], 0.0, atol=1e-12) or
            np.allclose(df[df.interaction_type == 'Enter Jaw R']['parent_x'], 0.0, atol=1e-12))

    assert (np.allclose(df[df.interaction_type == 'Exit Jaw']['parent_x'], 0.0, atol=1e-12) or 
            np.allclose(df[df.interaction_type == 'Exit Jaw']['ds'], line[tcp].length, atol=1e-12))

def test_impacts_crystal():
    coll = xc.EverestCrystal(length=0.002, material=xc.materials.SiliconCrystal, bending_angle=149e-6,
                        width=0.002, height=0.05, side='+', lattice='strip', jaw=0.001)

    num_part = int(50000)
    x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
    px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
    y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
    py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
    part = xp.Particles(x=x_init, px=px_init, y=y_init, py=py_init, delta=0, p0c=4e11)

    io_buffer = xt.new_io_buffer(capacity=int(2e7)) # 4-5 GB of memory
    coll.record_scatterings = True
    impacts_cry = xt.start_internal_logging(elements=[coll], io_buffer=io_buffer, capacity=io_buffer.capacity)
    df_crystal = impacts_cry.to_pandas()
    coll.track(part)
    part.sort(interleave_lost_particles=True)

    assert (np.allclose(df_crystal[df_crystal.interaction_type == 'Exit Jaw']['ds'], 0.0, atol=1e-12))