import numpy as np
import pandas as pd
import io


def load_SixTrack_colldb(filename, emitx, emity=None):
    return CollDB(emitx=emitx, emity=emity, sixtrack_file=filename)
    
class CollDB:
    def __init__(self, *, emitx, emity=None, sixtrack_file=None):
        self._emitx = emitx
        self._emity = emitx if emity is None else emity
        self._betagamma_rel = None
        if sixtrack_file is not None:
            self.load_SixTrack(sixtrack_file)
        else:
            self._colldb = None

    def show(self):
        return pd.DataFrame({
                'gap':             self.gap,
                's_center':        self.s_center,
                'angle':           self.angle,
                'material':        self.material,
                'offset':          self.offset,
                'tilt':            self.tilt,
                'stage':           self.stage,
                'active_length':   self.active_length,
                'collimator_type': self.collimator_type,
            }, index=self.name)
    
    @property
    def name(self):
        return self._colldb.index.values


    # TODO: - VALIDATION OF TYPES (e.g. material, stage)
    #       - IMPLEMENTATION OF TILT
    #       - CRYSTAL PROPERTIES: only valid if crystal == True (add mask to _set_property?)
    #       - icoll
    #       - show as __repr__
    
    # The CollDB class has the following fields (those marked
    # with an * are set automatically and cannot be overwritten):
    #   - name
    #   - gap
    #   - opening *
    #   - beam_size *
    #   - s_center *
    #   - angle
    #   - material
    #   - offset
    #   - tilt
    #   - stage
    #   - onesided
    #   - active_length
    #   - inactive_front
    #   - inactive_back
    #   - total_length *
    #   - collimator_type *
    #   - betx
    #   - bety
    #   - x
    #   - px
    #   - y
    #   - py
    #   - gamma_rel
    #   - emit
    #
    # Additionally, for crystals the following fields are added:
    #   - crystal
    #   - bend
    #   - xdim
    #   - ydim
    #   - miscut
    #   - thick
    
    @property
    def angle(self):
        return self._colldb['angle']

    @angle.setter
    def angle(self, angle):
        self._set_property('angle', angle)

    @property
    def material(self):
        return self._colldb['material']

    @material.setter
    def material(self, material):
        self._set_property('material', material)

    @property
    def offset(self):
        return self._colldb['offset']

    @offset.setter
    def offset(self, offset):
        self._set_property('offset', offset)

    @property
    def tilt(self):
        tilts = np.array([self._colldb.tilt_L.values,self._colldb.tilt_R.values])
        return pd.Series([ L if L == R else [L,R] for L, R in tilts.T ], index=self._colldb.index, dtype=object)

    @tilt.setter
    def tilt(self, tilts):
        self._set_property_LR('tilt', tilts)

    @property
    def stage(self):
        return self._colldb['stage']

    @stage.setter
    def stage(self, stage):
        self._set_property('stage', stage)

    @property
    def parking(self):
        return self._colldb['parking']

    @parking.setter
    def parking(self, parking):
        self._set_property('parking', parking)

#     @property
#     def crystal(self):
#         return self._colldb['crystal']

#     @crystal.setter
#     def crystal(self, crystal):
#         self._set_property('crystal', crystal)

#     @property
#     def bend(self):
#         return self._colldb['bend']

#     @bend.setter
#     def bend(self, bend):
#         self._set_property('bend', bend)

#     @property
#     def xdim(self):
#         return self._colldb['xdim']

#     @xdim.setter
#     def xdim(self, xdim):
#         self._set_property('xdim', xdim)

#     @property
#     def ydim(self):
#         return self._colldb['ydim']

#     @ydim.setter
#     def ydim(self, ydim):
#         self._set_property('ydim', ydim)

#     @property
#     def miscut(self):
#         return self._colldb['miscut']

#     @miscut.setter
#     def miscut(self, miscut):
#         self._set_property('miscut', miscut)

#     @property
#     def thick(self):
#         return self._colldb['thick']

#     @thick.setter
#     def thick(self, thick):
#         self._set_property('thick', thick)

    @property
    def s_center(self):
        return self._colldb['s_center']

    @property
    def collimator_type(self):
        return self._colldb['collimator_type']
    
    @property
    def active_length(self):
        return self._colldb['active_length']

    @active_length.setter
    def active_length(self, lengths):
        self._set_property('active_length', lengths)

    @property
    def inactive_front(self):
        return self._colldb['inactive_front']
    
    @inactive_front.setter
    def inactive_front(self, lengths):
        self._set_property('inactive_front', lengths)

    @property
    def inactive_back(self):
        return self._colldb['inactive_back']
    
    @inactive_back.setter
    def inactive_back(self, lengths):
        self._set_property('inactive_back', lengths)

    @property
    def total_length(self):
        return self._colldb['active_length'] +  self._colldb['inactive_front'] + self._colldb['inactive_back']
    
    @property
    def gap(self):
        gaps = np.array([self._colldb.gap_L.values,self._colldb.gap_R.values])
        return pd.Series([ L if L == R else [L,R] for L, R in gaps.T ], index=self._colldb.index, dtype=object)
    
    @gap.setter
    def gap(self, gaps):
        df = self._colldb
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
                mask_L = df.onesided.isin(['both','left'])
                mask_R = df.onesided.isin(['both','right'])
                df.loc[mask_L, 'gap_L'] = gaps[mask_L]
                df.loc[~mask_L, 'gap_L'] = None
                df.loc[mask_R, 'gap_R'] = gaps[mask_R]
                df.loc[~mask_R, 'gap_R'] = None
        
        # The variable gaps is a dictionary
        if isinstance(gaps, dict):
            correct_format = True
            for name, gap in gaps.items():
                if name not in self.name:
                    raise ValueError(f"Collimator {name} not found in CollDB!")
                side = df.onesided[name]
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
                df.loc[name, 'gap_L'] = gap_L
                df.loc[name, 'gap_R'] = gap_R

        if not correct_format:
            raise ValueError("Variable `gaps` needs to be a pandas Series, dict, numpy array, or list!")
        
        df.gap_L = df.gap_L.astype('object', copy=False)
        df.gap_R = df.gap_R.astype('object', copy=False)
            
        # Recompute openings
        self._compute_opening

    @property
    def opening(self):
        openings = np.array([self._colldb.opening_L.values,self._colldb.opening_R.values])
        return pd.Series([ L if L == R else [L,R] for L, R in openings.T ], index=self._colldb.index, dtype=object)
    
    @property
    def onesided(self):
        return self._colldb.onesided
    
    @onesided.setter
    def onesided(self, sides):
        self._set_property('onesided', sides)
        df = self._colldb
        df['onesided'] = [ 'both' if s==0 else s for s in df['onesided'] ]    
        df['onesided'] = [ 'left' if s==1 else s for s in df['onesided'] ]
        df['onesided'] = [ 'right' if s==2 else s for s in df['onesided'] ]
        self.gap = self.gap

    @property
    def gamma_rel(self):
        return np.sqrt(self._beta_gamma_rel**2+1)

    @gamma_rel.setter
    def gamma_rel(self, gamma_rel):
        self._beta_gamma_rel = np.sqrt(gamma_rel**2-1)
        self._compute_opening()

    @property
    def emit(self):
        return [self._emitx, self._emity]
    
    @emit.setter
    def emit(self, emit):
        if hasattr(emit, '__iter__'):
            if isinstance(emit, str):
                raise ValueError(f"The `emit` setting has to be a number!")
            elif len(emit) == 2:
                self._emitx = emit[0]
                self._emity = emit[1]
            elif len(val) == 1:
                self._emitx = emit[0]
                self._emity = emit[0]
            else:
                raise ValueError(f"The `emit` setting must have one or two values (for emitx and emity)!")
        else:
            self._emitx = emit
            self._emity = emit
        self._compute_opening()

    @property
    def betx(self):
        return self._colldb['betx']
    
    @betx.setter
    def betx(self, betx):
        self._set_property('betx', betx)
        self._compute_opening()

    @property
    def bety(self):
        return self._colldb['bety']
    
    @bety.setter
    def bety(self, bety):
        self._set_property('bety', bety)
        self._compute_opening()

    @property
    def x(self):
        return self._colldb['x']
    
    @x.setter
    def x(self, x):
        self._set_property('x', x)

    @property
    def px(self):
        return self._colldb['px']
    
    @px.setter
    def px(self, px):
        self._set_property('px', px)

    @property
    def y(self):
        return self._colldb['y']
    
    @y.setter
    def y(self, y):
        self._set_property('y', y)

    @property
    def py(self):
        return self._colldb['py']
    
    @py.setter
    def py(self, py):
        self._set_property('py', py)

    @property
    def beam_size(self):
        df = self._colldb
        return np.sqrt(
                (df['sigmax']*np.cos(np.float_(df['angle'].values)*np.pi/180))**2
                + (df['sigmay']*np.sin(np.float_(df['angle'].values)*np.pi/180))**2
            )

    def _compute_opening(self):
        df = self._colldb
        incomplete = np.any([ np.any([ x is None for x in df[opt] ]) for opt in ['betx', 'bety'] ])
        if self._beta_gamma_rel is None or incomplete:
            pass
        else:
            df['sigmax'] = np.sqrt(df['betx']*self._emitx/self._beta_gamma_rel)
            df['sigmay'] = np.sqrt(df['bety']*self._emity/self._beta_gamma_rel)
            df['opening_L'] = df['parking'] if df['gap_L'] is None else np.minimum(df['gap_L']*self.beam_size,df['parking'])
            df['opening_R'] = df['parking'] if df['gap_R'] is None else np.minimum(df['gap_R']*self.beam_size,df['parking'])



    # ---------------------------------------
    # ------ Property setter functions ------
    # ---------------------------------------

    def _set_property(self, prop, vals):
        df = self._colldb
        if isinstance(vals, dict):
            for name, val in vals.items():
                if name not in self.name:
                    raise ValueError(f"Collimator {name} not found in CollDB!")
                df.loc[name, prop] = val
        elif isinstance(vals, pd.Series) or isinstance(vals, list) or isinstance(vals, np.ndarray):
            if len(vals) != len(self.name):
                raise ValueError(f"The variable `{prop}` has a different length than the number of collimators in the CollDB. "
                                + "Use a dictionary instead.")
            df[prop] = vals
        else:
            raise ValueError(f"Variable `{prop}` needs to be a pandas Series, dict, numpy array, or list!")


    def _set_property_LR(self, prop, vals):
        df = self._colldb
        correct_format = False
        # The variable vals is a Series or a list
        if isinstance(vals, pd.Series) or isinstance(vals, list) or isinstance(vals, np.ndarray):
            correct_format = True
            if len(vals) != len(self.name):
                raise ValueError(f"The variable `{prop}` has a different length than the number of collimators in the CollDB. "
                                + "Use a dictionary instead.")
            # Some of the vals are list (e.g. two different values for both gaps): loop over vals as dict
            if any(hasattr(val, '__iter__') for val in vals):
                vals = dict(zip(self.name, vals))
            # All gaps are single values: use pandas-style assignment
            else:
                df[prop + "_L"] = vals
                df[prop + "_R"] = vals
        
        # The variable gaps is a dictionary
        if isinstance(vals, dict):
            correct_format = True
            for name, val in vals.items():
                if name not in self.name:
                    raise ValueError(f"Collimator {name} not found in CollDB!")
                if hasattr(val, '__iter__'):
                    if isinstance(val, str):
                        raise ValueError(f"The `{prop}` setting has to be a number!")
                    elif len(val) == 2:
                        val_L = val[0]
                        val_R = val[1]
                    elif len(val) == 1:
                        val_L = val[0]
                        val_R = val[0]
                    else:
                        raise ValueError(f"The `{prop}` setting must have one or two values (for the left and the right jaw)!")
                else:
                    val_L = val
                    val_R = val
                df.loc[name, prop + "_L"] = val_L
                df.loc[name, prop + "_R"] = val_R

        if not correct_format:
            raise ValueError("Variable `{prop}` needs to be a pandas Series, dict, numpy array, or list!")

        df[prop + "_L"] = df[prop + "_L"].astype('object', copy=False)
        df[prop + "_R"] = df[prop + "_R"].astype('object', copy=False)
                
            
    def _initialise_None(self):
        fields = {'s_center':None, 'gap_L': None, 'gap_R': None, 'opening_L': None, 'opening_R': None, 'angle': 0, 'offset': 0}
        fields.update({'parking': 0.025})
        fields.update({'onesided': 'both', 'material': None, 'stage': None, 'collimator_type': None, 'tilt_L': 0, 'tilt_R': 0})
        fields.update({'active_length': 0, 'inactive_front': 0, 'inactive_back': 0, 'sigmax': None, 'sigmay': None})
        fields.update({'crystal': None, 'bend': None, 'xdim': 0, 'ydim': 0, 'miscut': 0, 'thick': 0})
        fields.update({'betx': None, 'bety': None, 'x': None, 'px': None, 'y': None, 'py': None})
        for f, val in fields.items():
            if f not in self._colldb.columns:
                self._colldb[f] = val





    # -------------------------------
    # ------ Loading functions ------
    # -------------------------------

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
        
        sides = df['name'].apply(lambda s: onesided.get(s, 0))
        gaps = df.astype('object', copy=False)['opening'].apply(
            lambda s: None if float(family_settings.get(s, s)) >= 900 else float(family_settings.get(s, s))
        ).astype('object', copy=False)

        df['name'] = df['name'].str.lower() # Make the names lowercase for easy processing
        df.rename(columns={'length':'active_length'}, inplace=True)
        df = df.set_index('name')
        self._colldb = df.drop('opening', 1)
        
        self._initialise_None()
        print([ type(gap) for gap in gaps.values])
        self.gap = gaps.values
        self.onesided = sides.values

        # Check if collimator active
        # Check if gap is list (assymetric jaws)

        return