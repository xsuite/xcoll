import json
from pathlib import Path

import sixtrack_files as sf
import numpy as np
import xobjects as xo
import xpart as xp
import xcoll as xc

p0c = np.sqrt(6.8e12**2 - xp.PROTON_MASS_EV**2)

with open(Path(Path.cwd(), 'initial.json'), 'r') as fid:
    part_init = xp.Particles.from_dict(json.load(fid))

path = Path(Path.cwd(), 'SixTrack_B1')

with open(Path(path, 'collimators'), 'r') as fid:
    collimators = [ x.strip() for x in fid.readlines()]
for name in collimators:
    file = Path(path,'output','dump_mt_' + name + '_mken')
    dump_in = sf.sixtrack_dump2_to_particles(file, p0c=p0c)
    assert np.allclose(part_init.x,     dump_in.x, atol=1e-11, rtol=0)
    assert np.allclose(part_init.y,     dump_in.y, atol=1e-11, rtol=0)
    assert np.allclose(part_init.px,    dump_in.px, atol=1e-15, rtol=0)
    assert np.allclose(part_init.py,    dump_in.py, atol=1e-15, rtol=0)
    assert np.allclose(part_init.zeta,  dump_in.zeta, atol=1e-10, rtol=0)
    assert np.allclose(part_init.delta, dump_in.delta, atol=1e-13, rtol=0)
    file = Path(path,'output','dump_mt_' + name + '_mkex')
    dump_out = sf.sixtrack_dump2_to_particles(file, p0c=p0c)
    file = Path(Path.cwd(),'Ref',name + '.json')
    with open(file, 'w') as fid:
        json.dump(dump_out.to_dict(), fid, cls=xo.JEncoder)


path = Path(Path.cwd(), 'SixTrack_B2')

with open(Path(path, 'collimators'), 'r') as fid:
    collimators = [ x.strip() for x in fid.readlines()]
for name in collimators:
    file = Path(path,'output','dump_mt_' + name + '_mken')
    dump_in = sf.sixtrack_dump2_to_particles(file, p0c=p0c)
    assert np.allclose(part_init.x,     dump_in.x, atol=1e-11, rtol=0)
    assert np.allclose(part_init.y,     dump_in.y, atol=1e-11, rtol=0)
    assert np.allclose(part_init.px,    dump_in.px, atol=1e-15, rtol=0)
    assert np.allclose(part_init.py,    dump_in.py, atol=1e-15, rtol=0)
    assert np.allclose(part_init.zeta,  dump_in.zeta, atol=1e-10, rtol=0)
    assert np.allclose(part_init.delta, dump_in.delta, atol=1e-13, rtol=0)
    file = Path(path,'output','dump_mt_' + name + '_mkex')
    dump_out = sf.sixtrack_dump2_to_particles(file, p0c=p0c)
    file = Path(Path.cwd(),'Ref',name + '.json')
    with open(file, 'w') as fid:
        json.dump(dump_out.to_dict(), fid, cls=xo.
