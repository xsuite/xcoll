import numpy as np
import matplotlib.pyplot as plt
import xtrack as xt
import xpart as xp
import xobjects as xo
import xcoll as xc


env = xt.load(xc._pkg_root.parent / 'examples' / 'machines' / 'sps_q20_inj.json')
tt = env.sps.get_table()
tt = tt.rows[[not nn.startswith('drift_') for nn in tt.name]]

line = xt.load(xc._pkg_root.parent / 'examples' / 'machines' / 'sps_with_aperture_inj_q20.json')
env_aper = line.env
env_aper['sps'] = line
tt_aper = env_aper.sps.get_table()
tt_aper = tt_aper.rows[[et.startswith('Limit') for et in tt_aper.element_type]]

# Verify all apertures are named as xxxxxxxxx.[abcdefz]_aper
print(np.unique([aa[:-5].split('.')[-1] for aa in tt_aper.name]))
