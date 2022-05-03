from pathlib import Path
import json
import numpy as np

import xobjects as xo
import xcoll as xc

with open(Path(Path.cwd(), 'collimators_B1'), 'r') as fid:
    collimators = [ x.strip() for x in fid.readlines()]
for name in collimators:
    with open(Path(Path.cwd(), 'Collimators',name+'.json'), 'r') as fid:
        coll = xc.K2Collimator.from_dict(json.load(fid))
    coll.dx  *= 1e-3
    coll.dy  *= 1e-3
    coll.dpx *= 1e-3
    coll.dpy *= 1e-3
    coll.angle *= 180./np.pi
    with open(Path(Path.cwd(), 'Collimators',name+'.json'), 'w') as fid:
        json.dump(coll.to_dict(), fid, cls=xo.JEncoder)

with open(Path(Path.cwd(), 'collimators_B2'), 'r') as fid:
    collimators = [ x.strip() for x in fid.readlines()]
for name in collimators:
    with open(Path(Path.cwd(), 'Collimators',name+'.json'), 'r') as fid:
        coll = xc.K2Collimator.from_dict(json.load(fid))
    coll.dx  *= 1e-3
    coll.dy  *= 1e-3
    coll.dpx *= 1e-3
    coll.dpy *= 1e-3
    coll.angle *= 180./np.pi
    with open(Path(Path.cwd(), 'Collimators',name+'.json'), 'w') as fid:
        json.dump(coll.to_dict(), fid, cls=xo.JEncoder)
