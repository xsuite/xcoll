#!/usr/bin/env python

import json
import sys
from pathlib import Path
import numpy as np

from cpymad.madx import Madx

import xobjects as xo
import xtrack as xt
import xpart as xp


input_file  = Path(str(sys.argv[1]))
beam        = int(sys.argv[2])

output_file = input_file.parent / (input_file.stem + '.json')


mad = Madx()
mad.call(input_file.as_posix())
sequence = 'lhcb1' if beam==1 else 'lhcb2'

line = xt.Line.from_madx_sequence(mad.sequence[sequence], apply_madx_errors=True, install_apertures=True)
print(f"Imported {len(line.element_names)} elements.")
line.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, pdg_id=2212, gamma0=mad.sequence[sequence].beam.gamma)


# Elements to keep
collimators = [name for name in line.element_names
                    if (name.startswith('tc') or name.startswith('td'))
                    and not '_aper' in name and not name[-4:-2]=='mk' and not name[:4] == 'tcds'
                    and not name[:4] == 'tcdd' and not name[:5] == 'tclim' and not name[:3] == 'tca'
                    and not (name[-5]=='.' and name[-3]=='.') and not name[:5] == 'tcdqm'
              ]
# collimator_apertures = [f'{coll}_aper' + p for p in ['', '_patch'] for coll in collimators]
# ips = [f'ip{i+1}' for i in range(8)]


# Patch the aperture model by fixing missing apertures
df = line.check_aperture(needs_aperture=collimators)
missing_apertures = df.loc[df.has_aperture_problem, 'name'].values
def patch_aperture_with(missing, patch):
    if isinstance(missing, str) or not hasattr(missing, '__iter__'):
        missing = [missing]
    for nn in missing:
        if nn not in missing_apertures:
            print(f"No need to patch {nn}, aperture is present")
            continue
        if nn not in line.element_names:
            print(f"Element {nn} not found in line! Skipping aperture patching..")
            continue
        if isinstance(patch, str):
            if patch not in line.element_names:
                raise ValueError("Could not find patch aperture!")
            patch = line[patch].copy()
        line.insert_element(index=nn, element=patch,
                            name=nn+'_aper_patch')

if beam == 1:
    patch_aperture_with(['mo.28r3.b1', 'mo.32r3.b1'], 'mo.22r1.b1_mken_aper')
    patch_aperture_with(['mqwa.f5l7.b1..1', 'mqwa.f5l7.b1..2', 'mqwa.f5l7.b1..3',
                         'mqwa.f5l7.b1..4', 'mqwa.f5r7.b1..1', 'mqwa.f5r7.b1..2',
                         'mqwa.f5r7.b1..3', 'mqwa.f5r7.b1..4'
                        ], 'mqwa.e5l3.b1_mken_aper')
    patch_aperture_with(['tdisa.a4l2.b1', 'tdisb.a4l2.b1', 'tdisc.a4l2.b1'
                        ], xt.LimitRect(min_x=-0.043, max_x=0.043, min_y=-0.055, max_y=0.055))
    patch_aperture_with('tcld.a11r2.b1', xt.LimitEllipse(a=4e-2, b=4e-2))
    patch_aperture_with(['tcspm.b4l7.b1', 'tcspm.e5r7.b1', 'tcspm.6r7.b1'
                        ], xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-06))
    patch_aperture_with(['tcpch.a4l7.b1', 'tcpcv.a6l7.b1'
                        ], xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-06))

else:
    patch_aperture_with(['mo.32r3.b2', 'mo.28r3.b2'], 'mo.22l1.b2_mken_aper')
    patch_aperture_with(['mqwa.f5r7.b2..1', 'mqwa.f5r7.b2..2', 'mqwa.f5r7.b2..3',
                         'mqwa.f5r7.b2..4', 'mqwa.f5l7.b2..1', 'mqwa.f5l7.b2..2',
                         'mqwa.f5l7.b2..3', 'mqwa.f5l7.b2..4'
                        ], 'mqwa.e5r3.b2_mken_aper')
    patch_aperture_with(['tdisa.a4r8.b2', 'tdisb.a4r8.b2', 'tdisc.a4r8.b2'
                        ], xt.LimitRect(min_x=-0.043, max_x=0.043, min_y=-0.055, max_y=0.055))
    patch_aperture_with('tcld.a11l2.b2', xt.LimitEllipse(a=4e-2, b=4e-2))
    patch_aperture_with(['tcspm.d4r7.b2', 'tcspm.b4r7.b2', 'tcspm.e5l7.b2', 'tcspm.6l7.b2'
                        ], xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-06))
    patch_aperture_with(['tcpch.a5r7.b2', 'tcpcv.a6r7.b2'
                        ], xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-06))


# Do some simplifications on the line
# No, we don't want to do this by default
# line.remove_inactive_multipoles(inplace=True)
# line.remove_redundant_apertures(inplace=True, keep=collimator_apertures)
# line.remove_markers(inplace=True, keep=collimators+ips)
# line.remove_zero_length_drifts(inplace=True, keep=collimators)
# line.merge_consecutive_drifts(inplace=True, keep=collimators)
# print(f"Reduced to {len(line.element_names)} elements!")


print("Aperture check after patching:")
df_patched = line.check_aperture(needs_aperture=collimators)
assert not np.any(df_patched.has_aperture_problem)


# Save to json
with open(output_file, 'w') as fid:
    json.dump(line.to_dict(), fid, cls=xo.JEncoder, indent=True)


# Load from json to check that there are no loading errors
print("Reloading file to test json is not corrupted..")
with open(output_file, 'r') as fid:
    loaded_dct = json.load(fid)
newline = xt.Line.from_dict(loaded_dct)
assert xt._lines_equal(line, newline)
print("All done.")





