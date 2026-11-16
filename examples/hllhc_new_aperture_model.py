from pathlib import Path
import xtrack as xt
import xcoll as xc
import numpy as np

def sandwich_elements_by_aperture(line, *, need_aperture=[], tt=None):
    if tt is None:
        tt = line.get_table()
    if not line.aperture.is_built:
        line.aperture.build_model()
    aper_s = []
    aper_n = []
    insertions = []
    for nn, ss, ee, et in zip(tt.name, tt.s_start, tt.s_end, tt.element_type):
        if nn not in need_aperture:
            if et.startswith('Marker') or et.startswith('Limit') \
            or et.startswith('Drift') or et.startswith('Device'):
                continue
        if nn == '_end_point':
            continue
        aper_s.append([ss, ee])
        aper_n.append(nn)
    aper_s = np.array(aper_s)
    aper = line.aperture.get_limit_elements_at_s(aper_s.flatten())
    for name, (aper1, aper2) in zip(aper_n, aper_s):
        nn1 = f'{name}/aper_start'
        nn2 = f'{name}/aper_end'
        line.env.elements[nn1] = aper[aper1]
        line.env.elements[nn2] = aper[aper2]
        insertions.append(line.env.place(nn1, at=name+'@start'))
        insertions.append(line.env.place(nn2, at=name+'@end'))
    line.insert(insertions)


path = Path.cwd().parent / 'examples'
lhc = xt.load('https://acc-models.web.cern.ch/acc-models/lhc/hl19/xsuite/lhc.json')
line = lhc['b1']

lhc.aperture.from_json('https://acc-models.web.cern.ch/acc-models/lhc/hl19/aperture/aperture.json')

tt = line.get_table()
need_aperture = tt.name[tt.base_prototype == 'collimator']
sandwich_elements_by_aperture(line, need_aperture=need_aperture, tt=tt)

colldb = xc.CollimatorDatabase.from_yaml(path / 'colldbs' / 'hllhc_tight_15cm.yaml', beam=1)
colldb._collimator_dict = {
    kk.replace('tclia.4r2','tclia.4r2/b1'): vv
    for kk, vv in colldb._collimator_dict.items()
}
colldb.install_everest_collimators(verbose=True, line=line)
