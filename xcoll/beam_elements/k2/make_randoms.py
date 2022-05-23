try:
    import xcoll.beam_elements.pyk2 as pyk2
except ImportError:
    raise Exception("Error: Failed importing pyK2 (did you compile?). Cannot track.")

import numpy  as np
# import json
from pathlib import Path


seed = 6574654
pyk2.pyk2_init(seed)

num = int(1e7)
pseudo_rands = np.zeros(num, dtype=np.float64)
for i in range (num):
    pseudo_rands[i] = pyk2.pyk2_rand()

np.save(Path(Path.cwd(),'randoms.npy'),pseudo_rands)

# pseudo_rands_1 = np.zeros(num, dtype=np.float64)
# for i in range (num):
#     xran = np.array(0, dtype=np.float64)
#     pseudo_rands_1[i] = pyk2.pyk2_funlux()

# np.save(Path(Path.cwd(),'randoms1.npy'),pseudo_rands_1)

