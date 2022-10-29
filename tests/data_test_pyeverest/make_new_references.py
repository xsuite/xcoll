import json
from pathlib import Path
import numpy as np

import xobjects as xo
import xpart as xp
import xcoll as xc

collimators = [
    'tcl.4r1.b1', 'tcl.5r1.b1', 'tcl.6r1.b1', 'tctph.4l2.b1', 'tcsg.5l3.b1', 'tcsg.4r3.b1', 'tcla.b5r3.b1', 'tcla.6r3.b1', \
    'tcla.7r3.b1', 'tctph.4l5.b1', 'tcl.4r5.b1', 'tcl.5r5.b1', 'tcl.6r5.b1', 'tcsp.a4r6.b1', 'tclia.4r2', 'tcsg.a5r3.b1', \
    'tcsg.b5r3.b1', 'tdisa.a4l2.b1', 'tcld.a11r2.b1', 'tcspm.b4l7.b1', 'tcspm.e5r7.b1', 'tcl.4l1.b2', 'tcl.5l1.b2', \
    'tcl.6l1.b2', 'tctph.4r8.b2', 'tcspm.b4r7.b2', 'tcla.b6l7.b2', 'tcla.d6l7.b2', 'tcla.a7l7.b2', 'tcsp.a4l6.b2', \
    'tctph.4r5.b2', 'tcl.4l5.b2', 'tcl.5l5.b2', 'tcl.6l5.b2', 'tcsg.5r3.b2', 'tctpv.4r8.b2', 'tcsg.a6r7.b2', 'tcsg.6l7.b2', \
    'tdisb.a4r8.b2', 'tcld.a11l2.b2', 'tcsg.4l3.b2', 'tcsg.b5l3.b2', 'tcp.c6l7.b1', 'tcp.c6r7.b2', 'tcpcv.a6l7.b1', \
    'tcpch.a4l7.b1', 'tcpcv.a6r7.b2', 'tcpch.a5r7.b2'
]

path = Path.cwd()

def _reshuffle(part):
    # part.move(_context=xo.ContextCpu())
    if part.lost_particles_are_hidden:
        part.unhide_lost_particles()
    sort = np.argsort(part.particle_id)
    with part._bypass_linked_vars():
        for tt, nn in part._structure['per_particle_vars']:
            vv = getattr(part, nn)
            vv[:] = vv[sort]

def _make_collimator_ref(name):
    with open(Path(path, 'initial.json'), 'r') as fid:
        part = xp.Particles.from_dict(json.load(fid))
    with open(Path(path, 'Collimators', name+'.json'), 'r') as fid:
        colldict = json.load(fid)
    xc.scattering_routines.pyeverest.random.set_random_seed(6574654)
    if colldict['__class__'] == 'PyCollimator':
        coll = xc.beam_elements.PyCollimator.from_dict(colldict)
    elif colldict['__class__'] == 'PyCrystal':
        coll = xc.beam_elements.PyCrystal.from_dict(colldict)
    coll.track(part)
    _reshuffle(part)
    with open(Path(path, 'Ref',name+'.json'), 'w') as fid:
        json.dump(part.to_dict(), fid, cls=xo.JEncoder)

for coll in collimators:
    _make_collimator_ref(coll)


