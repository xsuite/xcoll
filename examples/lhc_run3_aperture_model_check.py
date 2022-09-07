import json
import numpy    as np
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack   as xt
import xpart    as xp
import xcoll    as xc



# Make a context and get a buffer
context = xo.ContextCpu()         # For CPU
# context = xo.ContextCupy()      # For CUDA GPUs
# context = xo.ContextPyopencl()  # For OpenCL GPUs
buffer = context.new_buffer()



# Load from json
with open('machines/lhc_run3_b1.json', 'r') as fid:
    loaded_dct = json.load(fid)
line = xt.Line.from_dict(loaded_dct)

line['acsca.d5l4.b1'].frequency = 400e6
line['acsca.c5l4.b1'].frequency = 400e6
line['acsca.b5l4.b1'].frequency = 400e6
line['acsca.a5l4.b1'].frequency = 400e6
line['acsca.a5r4.b1'].frequency = 400e6
line['acsca.b5r4.b1'].frequency = 400e6
line['acsca.c5r4.b1'].frequency = 400e6
line['acsca.d5r4.b1'].frequency = 400e6

elements = line.elements
s_elements = np.array(line.get_s_elements())
element_types = list(map(lambda e: e.__class__.__name__, elements))

import pandas as pd

elements_df = pd.DataFrame({
    'element_type': element_types,
    's': s_elements,
    'name': line.element_names
})

elements_df['is_aperture'] = elements_df.element_type.map(lambda s: s.startswith('Limit'))
elements_df['i_aperture_upstream'] = np.nan
elements_df['s_aperture_upstream'] = np.nan
elements_df['i_aperture_downstream'] = np.nan
elements_df['s_aperture_downstream'] = np.nan

num_elements = len(line.element_names)

i_prev_aperture = elements_df[elements_df['is_aperture']].index[0]
i_next_aperture = 0

# TODO: check no aperture before i_prev_aperture

for iee in range(i_prev_aperture, num_elements):

    print(f'{iee}\t', end='\r', flush=True)

    if elements_df.loc[iee, 'element_type'] == 'Drift':
        continue

    if elements_df.loc[iee, 'element_type'] == 'XYShift':
        continue

    if elements_df.loc[iee, 'element_type'] == 'SRotation':
        continue

    if elements_df.loc[iee, 'is_aperture']:
        i_prev_aperture = iee
        continue

    if i_next_aperture < iee:
        for ii in range(iee, num_elements):
            if elements_df.loc[ii, 'is_aperture']:
                i_next_aperture = ii
                break

    elements_df.at[iee, 'i_aperture_upstream'] = i_prev_aperture
    elements_df.at[iee, 'i_aperture_downstream'] = i_next_aperture

    elements_df.at[iee, 's_aperture_upstream'] = elements_df.loc[i_prev_aperture, 's']
    elements_df.at[iee, 's_aperture_downstream'] = elements_df.loc[i_next_aperture, 's']

elements_df['misses_aperture_upstream'] = ((elements_df['s_aperture_upstream'] != elements_df['s'])
    & ~(np.isnan(elements_df['i_aperture_upstream'])))



