# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
from pathlib import Path

import numpy as np
import xpart as xp
import xcoll as xc
import xobjects as xo


path = xc._pkg_root.parent / 'tests' / 'data_test_everest'


materials_b1 = {
  'BE':   'tcl.4r1.b1',
  'AL':   'tcl.5r1.b1',
  'CU':   'tcl.6r1.b1',
  'W':    'tctph.4l2.b1',
  'PB':   'tcsg.5l3.b1',
  'C':    'tcsg.4r3.b1',
  'C2':   'tcla.b5r3.b1',
  'Si':   'tcla.6r3.b1',
  'Ge':   'tcla.7r3.b1',
  'MoGR': 'tctph.4l5.b1',
  'CuCD': 'tcl.4r5.b1',
  'Mo':   'tcl.5r5.b1',
  'Glid': 'tcl.6r5.b1',
  'Iner': 'tcsp.a4r6.b1',
}

angles_b1 = {
  90.0:  'tclia.4r2',
  170.7: 'tcsg.a5r3.b1',
  10.8:  'tcsg.b5r3.b1',
}

lengths_b1 = {
  1.565: 'tdisa.a4l2.b1',
  0.600: 'tcld.a11r2.b1',
}


materials_b2 = {
  'BE':   'tcl.4l1.b2',
  'AL':   'tcl.5l1.b2',
  'CU':   'tcl.6l1.b2',
  'W':    'tctph.4r8.b2',
  'PB':   'tcspm.b4r7.b2',
  'C':    'tcla.b6l7.b2',
  'C2':   'tcla.d6l7.b2',
  'Si':   'tcla.a7l7.b2',
  'Ge':   'tcsp.a4l6.b2',
  'MoGR': 'tctph.4r5.b2',
  'CuCD': 'tcl.4l5.b2',
  'Mo':   'tcl.5l5.b2',
  'Glid': 'tcl.6l5.b2',
  'Iner': 'tcsg.5r3.b2',
}

angles_b2 = {
  90.0:  'tctpv.4r8.b2',
  141.1: 'tcsg.a6r7.b2',
  0.5:   'tcsg.6l7.b2',
}

lengths_b2 = {
  1.565: 'tdisb.a4r8.b2',
  0.600: 'tcld.a11l2.b2',
}


crystals_b1 = [
  'tcpcv.a6l7.b1',
  'tcpch.a4l7.b1'
]

crystals_b2 = [
  'tcpcv.a6r7.b2',
  'tcpch.a5r7.b2'
]


test_context = xo.ContextCpu()
# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
# def test_primaries(test_context):
def test_primaries():
    _track_collimator('tcp.c6l7.b1', _context=test_context)
    _track_collimator('tcp.c6r7.b2', _context=test_context)

# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
# def test_materials_b1(test_context):
def test_materials_b1():
    for key, name in materials_b1.items():
        _track_collimator(name, _context=test_context)

# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
# def test_materials_b2(test_context):
def test_materials_b2():
    for key, name in materials_b2.items():
        _track_collimator(name, _context=test_context)

# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
# def test_angles_b1(test_context):
def test_angles_b1():
    for key, name in angles_b1.items():
        _track_collimator(name, _context=test_context)

# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
# def test_angles_b2(test_context):
def test_angles_b2():
    for key, name in angles_b2.items():
        _track_collimator(name, _context=test_context)

# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
# def test_lengths_b1(test_context):
def test_lengths_b1():
    for key, name in lengths_b1.items():
        _track_collimator(name, _context=test_context, atolz=2e-11)

# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
# def test_lengths_b2(test_context):
def test_lengths_b2():
    for key, name in lengths_b2.items():
        _track_collimator(name, _context=test_context, atolz=2e-11)

# @for_all_test_contexts(
#     excluding=('ContextCupy', 'ContextPyopencl')  # Rutherford RNG not on GPU
# )
# def test_crystals(test_context):
def test_crystals():
    for name in crystals_b1 + crystals_b2:
        _track_collimator(name, _context=test_context)


def _track_collimator(name, atolx=3e-9, atoly=3e-9, atolpx=5e-9, atolpy=5e-9, atolz=1e-11, atold=2e-8, _context=None):
    print(f"Testing {name}")
    if _context is None:
        _context = xo.ContextCpu()
#     _context._cffi_verbose = True
#     _context._compile_kernels_info = False
    with open(Path(path, 'initial.json'), 'r') as fid:
        part = xp.Particles.from_dict(json.load(fid), _context=_context)
    with open(Path(path, 'Collimators', name+'.json'), 'r') as fid:
        colldict = json.load(fid)
    if colldict['__class__'] == 'EverestCollimator':
        coll = xc.EverestCollimator.from_dict(colldict, _context=_context)
    elif colldict['__class__'] == 'EverestCrystal':
        coll = xc.EverestCrystal.from_dict(colldict, _context=_context)
#   TODO: how to get compiled source, not save_source_as because 1) it doesnt work and 2) the line numbers are very wrong
#   I want the option for the test generated code (like 197660d59cf04c979191e2a334c2c7a0.c) to be NOT deleted afterwards
#   Also, I want to pass the following compiler flags: -Wall -Wextra for more info.
#     coll.compile_kernels(particles_class=xp.Particles, save_source_as='test.c')
    coll.track(part)
    part.sort(interleave_lost_particles=True)
    with open(Path(path, 'Ref',name+'.json'), 'r') as fid:
        part_ref = xp.Particles.from_dict(json.load(fid))
    part_ref.sort(interleave_lost_particles=True)
    assert np.array_equal(part.particle_id[part.state<1], part_ref.particle_id[part_ref.state<1])
    assert np.allclose(part.x[part.state>0],     part_ref.x[part_ref.state>0], atol=atolx, rtol=0)
    assert np.allclose(part.y[part.state>0],     part_ref.y[part_ref.state>0], atol=atoly, rtol=0)
    assert np.allclose(part.px[part.state>0],    part_ref.px[part_ref.state>0], atol=atolpx, rtol=0)
    assert np.allclose(part.py[part.state>0],    part_ref.py[part_ref.state>0], atol=atolpy, rtol=0)
    assert np.allclose(part.zeta[part.state>0],  part_ref.zeta[part_ref.state>0], atol=atolz, rtol=0)
    assert np.allclose(part.delta[part.state>0], part_ref.delta[part_ref.state>0], atol=atold, rtol=0)
    # for p, pref, pid in zip(part.x[part.state>0],     part_ref.x[part_ref.state>0], part.particle_id[part.state>0]):
    #     if not np.allclose(p, pref, atol=atolx, rtol=0):
    #         print(f"{pid}   x    : {abs(p-pref):.12}")
    # for p, pref, pid in zip(part.y[part.state>0],     part_ref.y[part_ref.state>0], part.particle_id[part.state>0]):
    #     if not np.allclose(p, pref, atol=atoly, rtol=0):
    #         print(f"{pid}   y    : {abs(p-pref):.12}")
    # for p, pref, pid in zip(part.px[part.state>0],    part_ref.px[part_ref.state>0], part.particle_id[part.state>0]):
    #     if not np.allclose(p, pref, atol=atolpx, rtol=0):
    #         print(f"{pid}   px   : {abs(p-pref):.12}")
    # for p, pref, pid in zip(part.py[part.state>0],    part_ref.py[part_ref.state>0], part.particle_id[part.state>0]):
    #     if not np.allclose(p, pref, atol=atolpy, rtol=0):
    #         print(f"{pid}   py   : {abs(p-pref):.12}")
    # for p, pref, pid in zip(part.zeta[part.state>0],  part_ref.zeta[part_ref.state>0], part.particle_id[part.state>0]):
    #     if not np.allclose(p, pref, atol=atolz, rtol=0):
    #         print(f"{pid}   zeta : {abs(p-pref):.12}")
    # for p, pref, pid in zip(part.delta[part.state>0], part_ref.delta[part_ref.state>0], part.particle_id[part.state>0]):
    #     if not np.allclose(p, pref, atol=atold, rtol=0):
    #         print(f"{pid}   delta: {abs(p-pref):.12}")
    # assert False

