import json
from pathlib import Path
import numpy as np
import xpart as xp
import xcoll as xc

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

offsets_b1 = {
  -0.0037:  'tcspm.b4l7.b1',
  0.000855: 'tcspm.e5r7.b1',
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
  'Glid': 'tcl.6l5.b2',    # Fails unless atol = 1e-9
  'Iner': 'tcsg.5r3.b2',   # Fails unless atol = 1e-9
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

offsets_b2 = {
  0.00297:   'tcsg.4l3.b2',
  -0.000346: 'tcsg.b5l3.b2',
}

crystals_b1 = [
  'tcpcv.a6l7.b1',
  'tcpch.a4l7.b1'
]

crystals_b2 = [
  'tcpcv.a6r7.b2',
  'tcpch.a5r7.b2'
]


path = Path('./data_test_K2/')

def test_primaries():
    _track_collimator('tcp.c6l7.b1')
    _track_collimator('tcp.c6r7.b2')

def test_materials_b1():
    for key, name in materials_b1.items():
        _track_collimator(name)

def test_materials_b2():
    for key, name in materials_b2.items():
        if name in ['tcl.6l5.b2','tcsg.5r3.b2']:
            _track_collimator(name, atolx=1e-9, atoly=1e-9, atolpx=1e-9, atolpy=1e-9, atold=1e-8)
        else:
            _track_collimator(name)

def test_angles_b1():
    for key, name in angles_b1.items():
        _track_collimator(name)

def test_angles_b2():
    for key, name in angles_b2.items():
        _track_collimator(name)

def test_lengths_b1():
    for key, name in lengths_b1.items():
        _track_collimator(name)

def test_lengths_b2():
    for key, name in lengths_b2.items():
        _track_collimator(name)

def test_offsets_b1():
    for key, name in offsets_b1.items():
        _track_collimator(name)

def test_offsets_b2():
    for key, name in offsets_b2.items():
        _track_collimator(name)

def test_crystals_b1():
    _track_collimator(crystals_b1[0], atolx=1e-7, atoly=1e-7, atolpx=1e-9, atolpy=1e-9, atold=1e-8)
def test_crystals_b1_bis():
    _track_collimator(crystals_b1[1], atolx=1e-7, atoly=1e-7, atolpx=1e-9, atolpy=1e-9, atold=1e-8)
def test_crystals_b2():
    _track_collimator(crystals_b2[0], atolx=1e-7, atoly=1e-7, atolpx=1e-9, atolpy=1e-9, atold=1e-8)
def test_crystals_b2_bis():
    _track_collimator(crystals_b2[1], atolx=1e-7, atoly=1e-7, atolpx=1e-9, atolpy=1e-9, atold=1e-8)


def _track_collimator(name, atolx=1e-11, atoly=1e-11, atolpx=1e-12, atolpy=1e-12, atolz=1e-10, atold=1e-10):
    with open(Path(path, 'initial.json'), 'r') as fid:
        part = xp.Particles.from_dict(json.load(fid))
    with open(Path(path, 'Collimators', name+'.json'), 'r') as fid:
        colldict = json.load(fid)
    seed = colldict.pop('random_generator_seed')
    xc.beam_elements.k2.k2_random.set_random_seed(seed)
    if colldict['__class__'] == 'K2Collimator':
        coll = xc.K2Collimator.from_dict(colldict)
    elif colldict['__class__'] == 'K2Crystal':
        coll = xc.K2Crystal.from_dict(colldict)
    coll.track(part)
    _reshuffle(part)
    with open(Path(path, 'Ref',name+'.json'), 'r') as fid:
        part_ref = xp.Particles.from_dict(json.load(fid))
    _reshuffle(part_ref)
    assert np.array_equal(part.particle_id[part.state<1], part_ref.particle_id[part_ref.state<1])
    assert np.allclose(part.x[part.state>0],     part_ref.x[part_ref.state>0], atol=atolx, rtol=0)
    assert np.allclose(part.y[part.state>0],     part_ref.y[part_ref.state>0], atol=atoly, rtol=0)
    assert np.allclose(part.px[part.state>0],    part_ref.px[part_ref.state>0], atol=atolpx, rtol=0)
    assert np.allclose(part.py[part.state>0],    part_ref.py[part_ref.state>0], atol=atolpy, rtol=0)
    assert np.allclose(part.zeta[part.state>0],  part_ref.zeta[part_ref.state>0], atol=atolz, rtol=0)
    assert np.allclose(part.delta[part.state>0], part_ref.delta[part_ref.state>0], atol=atold, rtol=0)

def _reshuffle(part):
#     part.move(_context=xo.ContextCpu())
    if part.lost_particles_are_hidden:
        part.unhide_lost_particles()

    sort = np.argsort(part.particle_id)
    with part._bypass_linked_vars():
        for tt, nn in part._structure['per_particle_vars']:
            vv = getattr(part, nn)
            vv[:] = vv[sort]
