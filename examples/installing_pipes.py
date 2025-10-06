import numpy as np
import matplotlib.pyplot as plt
import xtrack as xt
import xpart as xp
import xobjects as xo
import xcoll as xc


env = xt.load(xc._pkg_root.parent / 'examples' / 'machines' / 'sps_q20_inj.json')
tt = env.sps.get_table()
tt = tt.rows[[not nn.startswith('drift_') for nn in tt.name]]

# line = xt.load(xc._pkg_root.parent / 'examples' / 'machines' / 'sps_with_aperture_inj_q20.json')
# env_aper = line.env
# env_aper['sps'] = line
# tt_aper = env_aper.sps.get_table()
# tt_aper = tt_aper.rows[[et.startswith('Limit') for et in tt_aper.element_type]]

# !wget https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2025/aperture_ldb/APERTURE_EYETS%202024-2025.seq
with open('APERTURE_EYETS 2024-2025.seq', 'r') as f:
    mad_aper_content = f.read()
mad_aper_content = mad_aper_content.replace(
    "install, element = VEQF.10010.A, at= -0.095;",
   f"install, element = VEQF.10010.A, at= {line.get_length()-0.095};")

# Load aperture markers in a dummy sequence
from cpymad.madx import Madx
mad = Madx()
mad.input('''
SPS : SEQUENCE, refer = centre,    L = 6911.51818896;
a: marker, at = 20;
endsequence;
''')
mad.input(mad_aper_content)
mad.beam()
mad.use('SPS')
line_aper = xt.Line.from_madx_sequence(mad.sequence.SPS, install_apertures=True)
tt_aper = line_aper.get_table()
tt_aper = tt_aper.rows[[et.startswith('Limit') for et in tt_aper.element_type]]

# Verify all apertures are named as xxxxxxxxx.[abcdefz]_aper
print(np.unique([aa[:-5].split('.')[-1] for aa in tt_aper.name]))


class Pipe:
    def __init__(self, name, aperture_table):
        self.name = name
        self.table = aperture_table.rows[f"{name}.*"]

    @property
    def apertures(self):
        return [nn for _, nn in sorted(zip(self.positions, self.table.name))]

    @property
    def num_apertures(self):
        return len(self.apertures)

    @property
    def positions(self):
        return self.table.s

    @property
    def s_start(self):
        return min(self.positions)

    @property
    def s_end(self):
        return max(self.positions)

    @property
    def length(self):
        return self.s_end - self.s_start

    def __str__(self):
        return f"Pipe {self.name}: {self.num_apertures} apertures, length {self.length:.3f} m, from s={self.s_start:.3f} m to s={self.s_end:.3f} m"

    def __repr__(self):
        return self.__str__()

class ApertureModel:
    def __init__(self, aperture_table):
        self.table = aperture_table
        self.pipes = [Pipe(str(nn), aperture_table) for nn in np.unique([aa[:-7] for aa in aperture_table.name])]
        self.pipes = {pp.name: pp for pp in sorted(self.pipes, key=lambda pp: pp.s_start)}
        self._indices = {pp.name: ii for ii, pp in enumerate(self.pipes.values())}

    def get_overlaps(self, tol=1e-10):
        overlaps = []
        for ii in range(len(self.pipes)-1):
            pp1 = self[ii]
            pp2 = self[ii+1]
            if pp1.s_end > pp2.s_start:
                overlaps.append([pp1.name, pp2.name, float(pp1.s_end - pp2.s_start)])
        return overlaps

    def get_gaps(self, tol=1e-10):
        gaps = []
        for ii in range(len(self.pipes)-1):
            pp1 = self[ii]
            pp2 = self[ii+1]
            if pp1.s_end < pp2.s_start:
                gaps.append([pp1.name, pp2.name, float(pp2.s_start - pp1.s_end)])
        return gaps

    def index(self, name):
        try:
            return self._indices[name]
        except KeyError:
            raise ValueError(f"Pipe {name} not found")

    def __getitem__(self, key):
        if isinstance(key, (int, slice)):
            return list(self.pipes.values())[key]
        elif isinstance(key, str):
            return self.pipes[key]
        else:
            raise KeyError(f"Invalid key type: {type(key)}")

pipes = ApertureModel(tt_aper)
