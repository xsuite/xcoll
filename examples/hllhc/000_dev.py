import json
import io
import numpy as np
import pandas as pd

import xtrack as xt
import xpart as xp
import xcoll as xc

from temp_load_db import temp_load_colldb

print('Load file...')
with open('./HL_LHC_v1p5_clean_feb2022/HL_LHC_v1p5_line.json') as fid:
    dct = json.load(fid)
print('Build line...')
line = xt.Line.from_dict(dct)
#line = line.cycle(name_first_element='ip3')

# Attach reference particle (a proton a 7 TeV)
line.particle_ref = xp.Particles(mass0 = xp.PROTON_MASS_EV, p0c=7e12)



# Switch on RF (needed to twiss)
line['acsca.a5l4.b1'].voltage = 16e6
line['acsca.a5l4.b1'].frequency = 1e6

line0 = line.copy()

colldb = temp_load_colldb('HL_LHC_v1p5_clean_feb2022/CollDB_HL_relaxed_b1.data')

for kk in colldb.keys():
    assert kk in line.element_names

# Assumes marker in the the line is at the center of the active length
inputcolldf = pd.DataFrame.from_dict(colldb).transpose()\
                     .rename(columns={'length': 'active_length'})
inputcolldf['inactive_length_at_start'] = 1e-3
inputcolldf['inactive_length_at_end'] = 1e-3


parameters_to_be_extracted_from_twiss = (
    'x y px py betx bety alfx alfy gamx gamy dx dpx dy dpy mux muy'.split())
locations = ['at_center_active_part', 'at_start_active_part', 'at_start_element']

temp_dfs = []
for ll in locations:
    cols = pd.MultiIndex.from_tuples(
       [(ll, nn) for nn in ['s'] + parameters_to_be_extracted_from_twiss])
    newdf = pd.DataFrame(index=inputcolldf.index, columns=cols)
    temp_dfs.append(newdf)

colldf = temp_dfs[0]
for dd in temp_dfs[1:]:
    colldf = colldf.join(dd)

# Get collimator info from input
for cc in inputcolldf.columns[::-1]: #Try to get a nice ordering...
    colldf.insert(0, column=cc, value=inputcolldf[cc])

assert np.all(colldf['offset']==0), "Halfgap offset not implemented"

colldf['angle_rad'] = colldf['angle']
colldf['angle_deg'] = colldf['angle']*180/np.pi
colldf.drop('angle', axis=1, inplace=True, level=0)

colldf['length'] = (colldf['active_length']
                    + colldf['inactive_length_at_start']
                    + colldf['inactive_length_at_end'])
colldf['name'] = colldf.index

colldf['at_center_active_part', 's'] = line.get_s_position(colldf['name'].to_list())
colldf['at_start_active_part', 's'] = (
    colldf['at_center_active_part', 's'] - colldf['active_length']/2)
colldf['at_start_element', 's'] = (
    colldf['at_start_active_part', 's'] - colldf['inactive_length_at_start'])


for nn in colldf.index.values:
    print(nn)
    newcoll = xc.Collimator(
            inactive_length_at_start=colldf.loc[nn, 'inactive_length_at_start'],
            inactive_length_at_end=colldf.loc[nn, 'inactive_length_at_end'],
            active_length=colldf.loc[nn, 'active_length'],
            n_slices=10, angle=colldf.loc[nn, 'angle_deg'].values[0],
            a_min=-1, a_max=1,
            b_min=-1, b_max=1
            )
    print(newcoll.to_dict())
    s_insert = colldf['at_start_element', 's'][nn]
    line.insert_element(element=newcoll, name=nn, at_s=s_insert)

# Build tracker
tracker = xt.Tracker(line=line)

s_twiss = []
for ll in locations:
    s_twiss.extend(colldf[ll]['s'].to_list())

n_coll = len(colldf)
assert len(s_twiss) == len(locations) * n_coll

print('Start twiss')
tw = tracker.twiss(at_s=s_twiss)

# Extract optics info
for ill, ll in enumerate(locations):
    for nn in parameters_to_be_extracted_from_twiss:
        colldf[ll, nn] = tw[nn][ill*n_coll:(ill+1)*n_coll]

# Compute betatron beam sizes
nemitt_x_ref = 2.5e-6
nemitt_y_ref = 2.5e-6
beta0_gamma0 = (tracker.particle_ref._xobject.beta0[0]
                * tracker.particle_ref.gamma0[0])
for ll in locations:
    colldf[ll, 'sigmax'] = np.sqrt(colldf[ll, 'betx']*nemitt_x_ref/beta0_gamma0)
    colldf[ll, 'sigmay'] = np.sqrt(colldf[ll, 'bety']*nemitt_x_ref/beta0_gamma0)



# Compute halfgap
colldf['halfgap_m'] = colldf['nsigma'].values * np.sqrt(
      (colldf['at_center_active_part', 'sigmax']*np.cos(np.float_(colldf['angle_rad'].values)))**2
    + (colldf['at_center_active_part', 'sigmay']*np.sin(np.float_(colldf['angle_rad'].values)))**2)


# Configure collimators
for nn in colldf.index.values:
    line[nn].dx = colldf['at_center_active_part', 'x'][nn]
    line[nn].dy = colldf['at_center_active_part', 'y'][nn]
    line[nn].a_max = colldf['halfgap_m'][nn]
    line[nn].a_min = -colldf['halfgap_m'][nn]

# Machine aperture
n_sigmas = 30
n_part = 10000

x_norm = np.random.uniform(-n_sigmas, n_sigmas, n_part)
y_norm = np.random.uniform(-n_sigmas, n_sigmas, n_part)

part = xp.build_particles(tracker=tracker, x_norm=x_norm, y_norm=y_norm,
                          scale_with_transverse_norm_emitt=(2.5e-6, 2.5e-6),
                          at_element = 'ip3')

part0 = part.copy()

tracker.track(part, num_turns=5)


state_sorted = part.state.copy()
state_sorted[part.particle_id] = part.state
import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.plot(x_norm, y_norm, '.', color='red')
plt.plot(x_norm[state_sorted>0], y_norm[state_sorted>0], '.', color='green')
plt.axis('equal')
plt.show()
