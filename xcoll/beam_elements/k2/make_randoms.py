try:
    import xcoll.beam_elements.pyk2 as pyk2
except ImportError:
    raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")

import numpy  as np
# import json
from pathlib import Path
from xcoll.beam_elements.k2.materials import materials


seed = 6574654
pyk2.pyk2_init(seed)

num = int(1e7)
pseudo_rands = np.zeros(num, dtype=np.float64)
for i in range (num):
    pseudo_rands[i] = pyk2.pyk2_rand()

np.save(Path(Path.cwd(),'randoms.npy'),pseudo_rands)


for mat, val in materials.items():
    cgen = np.zeros(200, dtype=np.float64)
    zatom = val['zatom']
    emr = val['emr']
    hcut = val['hcut']
    pyk2.initialise_random(random_generator_seed=-1, cgen=cgen, zatom=zatom, emr=emr, hcut=hcut)
    np.save(Path(Path.cwd(),'cgen_' + mat + '.npy'),cgen)

