import json
import numpy as np

import xtrack as xt

print('Load file...')
with open('./HL_LHC_v1p5_clean_feb2022/HL_LHC_v1p5_line.json') as fid:
    dct = json.load(fid)

print('Build line...')
line = xt.Line.from_dict(dct)