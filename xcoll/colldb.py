import numpy as np
import pandas as pd
import io


def load_SixTrack_colldb(filename, emitx, emity=None):
    return CollDB(emitx=emitx, emity=emity, sixtrack_file=filename)
    
class CollDB:
    def __init__(self, *, emitx, emity=None, sixtrack_file=None):
        self.emitx = emitx
        self.emity = emitx if emity is None else emity
        self.gamma_rel = None
        if sixtrack_file is not None:
            self.load_SixTrack(sixtrack_file)
            self._colldb.insert(2,'s_center', None)
            self._colldb.insert(2,'opening', None)
            self._colldb['type'] = None
        else:
            self._colldb = None

    def show(self):
        return None
    
    @property
    def name(self):
        return self._colldb.index.values
    
    @property
    def gap(self):
        gaps = np.array([self._colldb.gap_L.values,self._colldb.gap_R.values])
        return pd.Series([ L if L == R else [L,R] for L, R in gaps.T ], index=self._colldb.index, dtype=object)
    
    @gap.setter
    def gap(self, gaps):
        correct_format = False
        # The variable gaps is a Series or a list
        if isinstance(gaps, pd.Series) or isinstance(gaps, list) or isinstance(gaps, np.ndarray):
            correct_format = True
            if len(gaps) != len(self.name):
                raise ValueError("The variable `gaps` has a different length than the number of collimators in the CollDB. "
                                + "Use a dictionary instead.")
            # Some of the gaps are list (e.g. two different values for both gaps): loop over gaps as dict
            if any(hasattr(gap, '__iter__') for gap in gaps):
                gaps = dict(zip(self.name, gaps))
            # All gaps are single values: use pandas-style assignment
            else:
                mask_L = self._colldb.onesided.isin(['both','left'])
                mask_R = self._colldb.onesided.isin(['both','right'])
                self._colldb.loc[mask_L, 'gap_L'] = gaps[mask_L]
                self._colldb.loc[~mask_L, 'gap_L'] = None
                self._colldb.loc[mask_R, 'gap_R'] = gaps[mask_R]
                self._colldb.loc[~mask_R, 'gap_R'] = None
        
        # The variable gaps is a dictionary
        if isinstance(gaps, dict):
            correct_format = True
            for name, gap in gaps.items():
                if name not in self.name:
                    raise ValueError(f"Collimator {name} not found in CollDB!")
                side = self._colldb.onesided[name]
                if hasattr(gap, '__iter__'):
                    if isinstance(gap, str):
                        raise ValueError("The gap setting has to be a number!")
                    elif len(gap) == 2:
                        gap_L = gap[0]
                        gap_R = gap[1]
                        if side != 'both':
                            if side == 'left' and gap_R is not None:
                                print(f"Warning: collimator {name} is left-sided but a finite right gap is specified. "
                                      + "Verify that this is what you want.")
                            elif side == 'right' and gap_L is not None:
                                print(f"Warning: collimator {name} is right-sided but a finite left gap is specified. "
                                      + "Verify that this is what you want.")
                    elif len(gap) == 1:
                        gap_L = gap[0] if side in ['both','left'] else None
                        gap_R = gap[0] if side in ['both','right'] else None
                    else:
                        raise ValueError("The gap setting must have one or two values (for the left and the right jaw)!")
                else:
                    gap_L = gap if side in ['both','left'] else None
                    gap_R = gap if side in ['both','right'] else None
                self._colldb.loc[name, 'gap_L'] = gap_L
                self._colldb.loc[name, 'gap_R'] = gap_R

        if not correct_format:
            raise ValueError("Variable `gaps` needs to be a pandas Series, dict, numpy array, or list!")
        
        self._colldb.gap_L = self._colldb.gap_L.astype('object', copy=False)
        self._colldb.gap_R = self._colldb.gap_R.astype('object', copy=False)
            
        # Recompute openings
        self._compute_opening


    @property
    def onesided(self):
        return self._colldb.onesided
    
    @onesided.setter
    def onesided(self, sides):
        df = self._colldb
        if isinstance(sides, dict):
            for name, side in sides.items():
                if name not in self.name:
                    raise ValueError(f"Collimator {name} not found in CollDB!")
                df.loc[name, 'onesided'] = side
        elif isinstance(sides, pd.Series) or isinstance(sides, list) or isinstance(sides, np.ndarray):
            if len(sides) != len(self.name):
                raise ValueError("The variable `sides` has a different length than the number of collimators in the CollDB. "
                                + "Use a dictionary instead.")
            df['onesided'] = sides
        else:
            raise ValueError("Variable `gaps` needs to be a pandas Series, dict, numpy array, or list!")
        df['onesided'] = [ 'both' if s==0 else s for s in df['onesided'] ]    
        df['onesided'] = [ 'left' if s==1 else s for s in df['onesided'] ]
        df['onesided'] = [ 'right' if s==2 else s for s in df['onesided'] ]
        self.gap = self.gap


    @property
    def sigma(self):
        sigma = np.array([self._colldb.sigmax.values,self._colldb.sigmay.values])
        return pd.Series(sigma.T, index=self._colldb.index, dtype=object) 

    @property
    def total_length(self):
        return self._colldb['active_length'] +  self._colldb['inactive_front'] + self._colldb['inactive_back']
    
    
    def _compute_opening(self):
        incomplete = np.any([ np.any([ x is None for x in colldb[opt] ]) for opt in ['betx', 'bety'] ])
        if self.gamma_rel is None or incomplete:
            pass
        else:
            self._colldb['sigmax'] = np.sqrt(self._colldb['betx']*self.emitx/self.gamma_rel)
            self._colldb['sigmay'] = np.sqrt(self._colldb['bety']*self.emity/self.gamma_rel)

            
            
    def _initialise_None(self):
        fields = {'gap_L': None, 'gap_R': None, 'opening_L': None, 'opening_R': None, 'active_length': 0}
        fields += {'angle': 0, 'onesided': 0, 'material': None, 'stage': None, 'offset': 0, 'tilt_L': 0, 'tilt_R': 0}
        fields += {'inactive_front': 0, 'inactive_back': 0}
        fields += {'crystal': None, 'bend': None, 'xdim': 0, 'ydim': 0, 'miscut': 0, 'thick': 0}
        fields += {'betx': None, 'bety': None, 'x': None, 'px': None, 'y': None, 'py': None, 'sigma_x': None, 'sigma_y': None}
        for f, val in fields.items():
            if f not in self._colldb.columns:
                self._colldb[f] = val
            
    def load_SixTrack(self,filename):
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

        df = df[['name', 'opening', 'length', 'angle', 'material', 'offset']]
        df.insert(5,'stage', df['opening'].apply(lambda s: family_types.get(s, 'UNKNOWN')))
        
        
        
        df['onesided'] = df['name'].apply(lambda s: onesided.get(s, 0))
        df['onesided'] = [ 'left' if s==1 else s for s in df['onesided'] ]
        df['onesided'] = [ 'right' if s==2 else s for s in df['onesided'] ]
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
        
        self._colldb = df.drop('opening', 1)
        self.gap = df['opening'].apply(lambda s: float(family_settings.get(s, s)))

        # Check if collimator active
        # Check if gap is list (assymetric jaws)

        return