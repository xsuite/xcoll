# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
import os
from pathlib import Path
import numpy as np

def rewrite_jaws(file):
    with open(file,'r') as fid:
        data = json.load(fid)
    s = (data['reference_center'][0]*np.cos(np.deg2rad(data['angle'])) + 
        data['reference_center'][1]*np.sin(np.deg2rad(data['angle']))) 

    if isinstance(data['jaw'], float): 
        data['jaw_L'] = s + data['jaw']
        data['jaw_R'] = s - data['jaw']
    else:
        data['jaw_L'] = s + data['jaw'][0]  
        data['jaw_R'] = s + data['jaw'][1]

    del data['reference_center']
    del data['jaw']

    with open(file, 'w') as fid:
        json.dump(data, fid, indent=1)

# Apply function to all the files in the directory
dir = 'Collimators'
files = os.listdir(dir)
for file in files:
    if file.endswith('.json'):
        rewrite_jaws(Path(dir) / file)
        print(f'{file} rewritten')
    else:
        print(f'{file} not rewritten')