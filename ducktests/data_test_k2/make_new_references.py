import json
import sys
from pathlib import Path
import numpy as np

import xobjects as xo
import xpart as xp
sys.path.append(str(Path.cwd().resolve().parents[1]))
import duckcoll as dc

dc.scattering_routines.k2.K2Engine(_capacity=50000, random_generator_seed=6574654)

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


def _make_collimator_ref(name):
    # Initialise engine
    dc.scattering_routines.k2.K2Engine.reset()
    # Load initial particles
    with open(Path(path, 'initial.json'), 'r') as fid:
        part = xp.Particles.from_dict(json.load(fid))
    # Initialise collimator
    with open(Path(path, 'Collimators', name+'.json'), 'r') as fid:
        colldict = json.load(fid)
    if colldict['__class__'] == 'K2Collimator':
        coll = dc.K2Collimator.from_dict(colldict)
    elif colldict['__class__'] == 'K2Crystal':
        coll = dc.K2Crystal.from_dict(colldict)
    # Track
    coll.track(part)
    part.sort(interleave_lost_particles=True)
    # Store reference
    with open(Path(path, 'Ref',name+'.json'), 'w') as fid:
        json.dump(part.to_dict(), fid, cls=xo.JEncoder)

for coll in collimators:
    _make_collimator_ref(coll)


