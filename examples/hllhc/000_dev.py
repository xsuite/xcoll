import json
import io
import numpy as np
import pandas as pd

import xtrack as xt

from temp_load_db import temp_load_colldb

print('Load file...')
with open('./HL_LHC_v1p5_clean_feb2022/HL_LHC_v1p5_line.json') as fid:
    dct = json.load(fid)

print('Build line...')
line = xt.Line.from_dict(dct)

colldb = temp_load_colldb('HL_LHC_v1p5_clean_feb2022/CollDB_HL_relaxed_b1.data')

for kk in colldb.keys():
    assert kk in line.element_names

# Assumes marker in the the line is at the center of the active length
inputcolldf = pd.DataFrame.from_dict(colldb).transpose()\
                     .rename(columns={'length': 'active_length'})
inputcolldf['inactive_length_at_start'] = 1e-2
inputcolldf['inactive_length_at_end'] = 1e-2



temp_dfs = []
for location in ['at_center_active_part',
                 'at_start_active_part', 'at_start_element']:
    cols = pd.MultiIndex.from_tuples(
        [(location, nn) for nn in
            's x y px py betx bety alfx alfy gamx gamy dx dpx dy dpy'.split()])
    newdf = pd.DataFrame(index=inputcolldf.index, columns=cols)
    temp_dfs.append(newdf)

colldf = temp_dfs[0]
for dd in temp_dfs[1:]:
    colldf = colldf.join(dd)

for cc in inputcolldf.columns[::-1]: #Try to get a nice ordering...
    colldf.insert(0, column=cc, value=inputcolldf[cc])

colldf['length'] = (colldf['active_length']
                    + colldf['inactive_length_at_start']
                    + colldf['inactive_length_at_end'])
colldf['name'] = colldf.index

colldf['at_center_active_part', 's'] = line.get_s_position(colldf['name'].to_list())
colldf['at_start_active_part', 's'] = (
    colldf['at_center_active_part', 's'] - colldf['active_length']/2)
colldf['at_start_element', 's'] = (
    colldf['at_start_active_part', 's'] - colldf['inactive_length_at_start'])


