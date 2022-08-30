from pathlib import Path
import json
import numpy as np

import xobjects as xo
import xcoll as xc


print("\n\nCorrecting collimator JSONs.\nThis should be ran only once, to avoid double corrections!\n\n")

for beam in 'B1', 'B2':
    with open(Path(Path.cwd(), 'collimators_' + beam), 'r') as fid:
        collimators = [ x.strip() for x in fid.readlines()]
    for name in collimators:
        with open(Path(Path.cwd(), 'Collimators',name+'.json'), 'r') as fid:
            colldict = json.load(fid)
            if colldict['__class__'] == 'K2Collimator':
                coll = xc.K2Collimator.from_dict(colldict)
            elif colldict['__class__'] == 'K2Crystal':
                coll = xc.K2Crystal.from_dict(colldict)
            else:
                raise ValueException(f"Class {colldict['__class__']} not understood!")
        coll.dx  *= 1e-3
        coll.dy  *= 1e-3
        coll.dpx *= 1e-3
        coll.dpy *= 1e-3
        coll.angle *= 180./np.pi
        with open(Path(Path.cwd(), 'Collimators',name+'.json'), 'w') as fid:
            json.dump(coll.to_dict(), fid, cls=xo.JEncoder)
