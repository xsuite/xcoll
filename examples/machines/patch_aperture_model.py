import json

import xtrack as xt
import xobjects as xo

# Load from json
with open('lhc_run3_b1_from_mask.json', 'r') as fid:
    loaded_dct = json.load(fid)
line = xt.Line.from_dict(loaded_dct)

print("\nAperture check before patching:")
df_imported = line.check_aperture()


# Patch aperture model
aperture_for_octupoles = line['mo.22r1.b1_mken_aper'].copy()

for nn in ['mo.28r3.b1', 'mo.32r3.b1']:
    line.insert_element(index=nn, element=aperture_for_octupoles,
                        name=nn+'_aper_patch')

aperture_for_mqwa = line['mqwa.e5l3.b1_mken_aper'].copy()
for nn in ['mqwa.f5l7.b1..1', 'mqwa.f5l7.b1..2', 'mqwa.f5l7.b1..3',
           'mqwa.f5l7.b1..4', 'mqwa.f5r7.b1..1', 'mqwa.f5r7.b1..2',
           'mqwa.f5r7.b1..3', 'mqwa.f5r7.b1..4']:
    line.insert_element(index=nn, element=aperture_for_mqwa,
                        name=nn+'_aper_patch')

aperture_for_tdis = line['tdis.4l2.a.b1_aper'].copy()
for nn in ['tdisa.a4l2.b1', 'tdisb.a4l2.b1', 'tdisc.a4l2.b1']:
    line.insert_element(index=nn, element=aperture_for_tdis,
                        name=nn+'_aper')

aperture_for_tcld = xt.LimitEllipse(a=4e-2, b=4e-2)
for nn in ['tcld.a11r2.b1']:
    line.insert_element(index=nn, element=aperture_for_tdis,
                            name=nn+'_aper')

aperture_for_tcspm = line['tcspm.6r7.a.b1_aper'].copy()
for nn in ['tcspm.b4l7.b1', 'tcspm.e5r7.b1', 'tcspm.6r7.b1']:
    line.insert_element(index=nn, element=aperture_for_tcspm,
                        name=nn+'_aper')

print("\nAperture check after patching:")
df_patched = line.check_aperture()

# Save to json
with open('lhc_run3_b1.json', 'w') as fid:
    json.dump(line.to_dict(), fid, cls=xo.JEncoder, indent=True)