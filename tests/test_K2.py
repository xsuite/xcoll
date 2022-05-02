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
  'CUCD': 'tcl.4r5.b1',
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
  3.000: 'tcdqa.a4r6.b1',
}

offsets_b1 = {
  -0.0037:  'tcspm.b4l7.b1',
  0.000855: 'tcspm.e5r7.b1',
}

path = Path('./data_test_K2/')
    
def _track_collimator(name):
    with open(Path(path, 'initial.json'), 'r') as fid:
        part = xp.Particles.from_dict(json.load(fid))
    with open(Path(path, 'Collimators', name+'.json'), 'r') as fid:
        colldict = json.load(fid)
    coll = xc.K2Collimator.from_dict(colldict)
    coll.track(part)
    part.reshuffle()
    with open(Path(path, 'Ref',name+'.json'), 'r') as fid:
        part_ref = xp.Particles.from_dict(json.load(fid))
    part_ref.reshuffle()
    assert np.array_equal(part_ref.particle_id[part_ref.state<1], test.particle_id[test.state<1])
    assert np.allclose(part.x,     part_ref.x, atol=1e-11, rtol=0)
    assert np.allclose(part.y,     part_ref.y, atol=1e-11, rtol=0)
    assert np.allclose(part.px,    part_ref.px, atol=1e-15, rtol=0)
    assert np.allclose(part.py,    part_ref.py, atol=1e-15, rtol=0)
    assert np.allclose(part.zeta,  part_ref.zeta, atol=1e-10, rtol=0)
    assert np.allclose(part.delta, part_ref.delta, atol=1e-13, rtol=0)
        
        
        
        
        