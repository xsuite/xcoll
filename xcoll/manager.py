import io

import numpy as np
import pandas as pd
from .beam_elements import Collimator
from pyk2 import K2Collimator

class CollimatorManager:
    def __init__(self, *, line, coll_db_txt_file, emitx, emity):

        colldb = _load_colldb(coll_db_txt_file)
        colldb.insert(1,'s_center', None)
        colldb.insert(1,'opening', None)
        colldb['type'] = None
        
        if not np.all(colldb['offset']==0):
            raise NotImplementedError("Halfgap offset not implemented")

        self.colldb = colldb
        self.line = line
        self.emitx = emitx
        self.emity = emity

    @property
    def collimator_names(self):
        return list(self.colldb.index.values)


    def install_black_absorbers(self, names=None):
        line = self.line
        colldb = self.colldb
        if names is None:
            names = colldb.index.values
        mask = colldb.index.isin(names)
        for name in names:
            if name not in line.element_names:
                raise Exception(f"Collimator {name} not found in line!")
        if line.tracker is not None:
            raise Exception("Tracker already built! Please install collimators before building tracker!")

        colldb.loc[mask,'s_center'] = line.get_s_position(names)
        colldb.loc[mask,'type'] = 'BlackAbsorber'

        for name in names:
            print(f"Installing {name}")
            thiscoll = colldb.loc[name]
            newcoll = Collimator(
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    jaw_R=-1, jaw_L=1,
                    jaw_D=-1, jaw_U=1
                    )
            s_install = thiscoll['s_center'] - thiscoll['active_length']/2 - thiscoll['inactive_front']
            line.insert_element(element=newcoll, name=name, at_s=s_install)


    def install_k2_collimators(self, names=None):
        line = self.line
        colldb = self.colldb
        if names is None:
            names = colldb.index.values
        mask = colldb.index.isin(names)
        for name in names:
            if name not in line.element_names:
                raise Exception(f"Collimator {name} not found in line!")
        if line.tracker is not None:
            raise Exception("Tracker already built! Please install collimators before building tracker!")

        colldb.loc[mask,'s_center'] = line.get_s_position(names)
        colldb.loc[mask,'type'] = 'K2Collimator'

        for name in names:
            print(f"Installing {name}")
            thiscoll = colldb.loc[name]
            newcoll = K2Collimator(
                    k2_engine,
                    length,
                    rotation,
                    icoll,
                    aperture,
                    onesided,
                    dx,
                    dy,
                    dpx,
                    dpy,
                    inactive_front=thiscoll['inactive_front'],
                    inactive_back=thiscoll['inactive_back'],
                    active_length=thiscoll['active_length'],
                    angle=thiscoll['angle'],
                    jaw_R=-1, jaw_L=1,
                    jaw_D=-1, jaw_U=1
                    )
            s_install = thiscoll['s_center'] - thiscoll['active_length']/2 - thiscoll['inactive_front']
            line.insert_element(element=newcoll, name=name, at_s=s_install)


    def _compute_optics(self):
        line = self.line
        if line is None or line.tracker is None:
            raise Exception("Please build tracker before calling this method!")
        
        colldb = self.colldb
        colldb['betx'] = None
        colldb['bety'] = None
        colldb['x'] = None
        colldb['px'] = None
        colldb['y'] = None
        colldb['py'] = None
        colldb['sigmax'] = None
        colldb['sigmay'] = None
        
        tracker = line.tracker
        tw = tracker.twiss(at_s=colldb['s_center'])
        colldb['betx'] = tw['betx']
        colldb['bety'] = tw['bety']
        colldb['x'] = tw['x']
        colldb['px'] = tw['px']
        colldb['y'] = tw['y']
        colldb['py'] = tw['py']
        beta0_gamma0 = tracker.particle_ref._xobject.beta0[0] * tracker.particle_ref._xobject.gamma0[0]
        colldb['sigmax'] = np.sqrt(colldb['betx']*self.emitx/beta0_gamma0)
        colldb['sigmay'] = np.sqrt(colldb['bety']*self.emity/beta0_gamma0)


    def set_openings(self):
        line = self.line
        if line is None or line.tracker is None:
            raise Exception("Please build tracker before calling this method!")
        colldb = self.colldb
        if any([ x is None for x in colldb.type ]):
            raise ValueError("Some collimators have not yet been installed. "
                             + "Please install all collimators before setting the openings.")
        # Compute halfgap
        self._compute_optics()
        colldb['opening'] = colldb['gap'].values * np.sqrt(
            (colldb['sigmax']*np.cos(np.float_(colldb['angle'].values)*np.pi/180))**2
            + (colldb['sigmay']*np.sin(np.float_(colldb['angle'].values)*np.pi/180))**2
        )

        # Configure collimators
        for name in colldb.index:
            if isinstance(line[name], Collimator):
                line[name].dx = colldb['x'][name]
                line[name].dy = colldb['y'][name]
                line[name].jaw_R = -colldb['opening'][name]
                line[name].jaw_L = colldb['opening'][name]
            elif isinstance(line[name], K2Collimator):
                pass
            else:
                raise ValueError(f"Missing implementation for element type of collimator {name}!")


def _load_colldb(filename):
    with open(filename, 'r') as infile:
        coll_data_string = ''
        family_settings = {}
        family_types = {}
        onesided = {}

        for l_no, line in enumerate(infile):
            if line.startswith('#'):
                continue # Comment

            sline = line.split()
            if len(sline) < 6:
                if sline[0].lower() == 'nsig_fam':
                    family_settings[sline[1]] = float(sline[2])
                    family_types[sline[1]] = sline[3]
                elif sline[0].lower() == 'onesided':
                    onesided[sline[1]] = int(sline[2])
                elif sline[0].lower() == 'settings':
                    pass # Acknowledge and ignore this line
                else:
                    print(f"Unknown setting {line}")
            else:
                coll_data_string += line

    names = ['name', 'opening', 'material', 'length', 'angle', 'offset']

    df = pd.read_csv(io.StringIO(coll_data_string), delim_whitespace=True,
                    index_col=False, names=names)

    #df['angle'] = df['angle']
    df = df[['name', 'opening', 'length', 'angle', 'material', 'offset']]
    df.insert(5,'stage', df['opening'].apply(lambda s: family_types.get(s, 'UNKNOWN')))
    df.insert(0,'gap', df['opening'].apply(lambda s: float(family_settings.get(s, s))))
    df['onesided'] = df['name'].apply(lambda s: onesided.get(s, 0))
    df['tilt_left'] = 0
    df['tilt_right'] = 0
    df['inactive_front'] = 0
    df['inactive_back'] = 0
    df['total_length'] = df['length']
    df['crystal'] = None
    df['bend'] = None
    df['xdim'] = None
    df['ydim'] = None
    df['miscut'] = None
    df['thick'] = None

    df['name'] = df['name'].str.lower() # Make the names lowercase for easy processing
    df.rename(columns={'length':'active_length'}, inplace=True)
    df = df.set_index('name')

    return df.drop('opening', 1)
