import numpy as np
import pandas as pd
import io

def load_SixTrack_colldb(filename, *, emit):
    return CollimatorDatabase(emit=emit, sixtrack_file=filename)

# in colldb:
#      (gap [sigma] + offset [m] + tilt [deg]) or physical_opening [m]
#      angle
#      length (all three)
#      onesided     TODO: need better name
#      parking ??
#      material
#      stage
#      family
#      crystal

# in element:
#      jaw
#      ref
#      angle
#      length (all three)
#      onesided ??
#      material
#      is_active

# in neither:
#      align_to
#      s_center
#      collimator_type

# if physical_opening is used, other variables like tilt, opening, offset, .. are ignored
_coll_properties = {'active_length': 0,
                    'inactive_front': 0,
                    'inactive_back': 0,
                    'gap_L': None,
                    'gap_R': None,
                    'onesided': 'both',
                    'angle_L': 0,
                    'angle_R': 0,
                    'offset': 0,
                    'tilt_L': 0,
                    'tilt_R': 0,
                    'align_to': None, 
                    's_center': None,
                    'parking': 1,
                    'material': None,
                    'jaw_LU': None,
                    'jaw_RU': None,
                    'jaw_LD': None,
                    'jaw_RD': None,
                    'ref_xU': 0,
                    'ref_yU': 0,
                    'ref_xD': 0,
                    'ref_yD': 0,
                    'stage': None,
                    'family': None,
                    'collimator_type': None,
                    'is_active': True,
                    'crystal': False
                   }
_properties_no_setter = ['jaw_LU', 'jaw_RU', 'jaw_LD', 'jaw_RD', 'gap_L', 'gap_R',
                         'ref_xU', 'ref_yU', 'ref_xD', 'ref_yD', 'angle_L', 'angle_R',
                         'tilt_L', 'tilt_R', ]
_properties_in_element = ['jaw_LU', 'jaw_RU', 'jaw_LD', 'jaw_RD', 'angle_L', 'angle_R', 
                          'ref_xU', 'ref_yU', 'ref_xD', 'ref_yD', 'active_length',
                          'inactive_front', 'inactive_back']
_add_to_dict = ['angle', 'tilt', 'opening', 'physical_opening', 'reference_center']

_crystal_properties = {
                    'bend': None,
                    'xdim': 0,
                    'ydim': 0,
                    'miscut': 0,
                    'thick': 0
                }
_optics_vals = ['x', 'px', 'y', 'py', 'betx', 'bety', 'alfx', 'alfy', 'dx', 'dy']


# This creates a view on the settings of one collimator, kept in sync with the main database
class CollimatorSettings:

    def __init__(self, name, colldb=None, optics=None, element=None):
        if colldb is None:
            colldb = pd.DataFrame({**_coll_properties, **_crystal_properties}, index=[name])
        self._colldb = colldb
        if name not in colldb.index:
            raise ValueError(f"Collimator {name} not found in database!")
        self._optics  = optics
        self._element = element
        self._name = name
        self._jaws_manually_set = False

        # Automatically assign all properties from _collimator_properties
        # Now create @property's for each of them:
        for pp in set(_coll_properties.keys()).union(set(_crystal_properties.keys())):
            if not pp in self._colldb.columns:
                raise ValueError(f"Error in settings for {self.name}: "
                               + f"Could not find property {pp} in database column header!")
            if pp in _properties_no_setter:
                setattr(self.__class__, pp, property(_prop_fget(self, pp)))
            else:
                setattr(self.__class__, pp, property(_prop_fget(self, pp), _prop_fset(self, pp)))

    @property
    def name(self):
        return self._name

    def to_dict(self):
        props = list(_coll_properties.keys()) + list(_crystal_properties.keys()) + _add_to_dict
        all_properties = {pp: getattr(self,pp) for pp in props if pp not in _properties_no_setter}
        return {'name': self.name, **all_properties}

    # Some easy accessors to the LR / LRUD properties:
    # -----------------------------------------------

    @property
    def angle(self):
        return _get_LR(self, 'angle')

    @angle.setter
    def angle(self, val):
        _set_LR(self, 'angle', val)
        self._compute_jaws()

    @property
    def tilt(self):
        return _get_LR(self, 'tilt')

    @tilt.setter
    def tilt(self, val):
        _set_LR(self, 'tilt', val)
        self._compute_jaws()

    @property
    def opening(self):
        if self._jaws_manually_set:
            # TODO: update gap using beam size as from jaws
            pass
        else:
            return _get_LR(self, 'gap')

    @opening.setter
    def opening(self, val):
        _set_LR(self, 'gap', val, neg=True)
        self._jaws_manually_set = False
        self._compute_jaws()

    @property
    def physical_opening(self):
        return _get_LRUD(self, 'jaw', neg=True)

    @physical_opening.setter
    def physical_opening(self, val):
        _set_LRUD(self, 'jaw', val, neg=True)
        self._jaws_manually_set = True
        self._compute_jaws()

    @property
    def reference_center(self):
        return _get_LRUD(self, 'ref', name_LU='_xU', name_RU='_yU', name_LD='_xD', name_RD='_yD')

    @property
    def beam_size(self):
        if not self._optics_is_ready:
            return None
        else:
            # TODO: calculate beam size from optic
            # [[LU,RU], [LC, RC], [LD,RD]]
            pass


    @property
    def _optics_is_ready(self):
        if self.align_to is None or self._optics is None:
            return False
        # TODO: when True

    def _compute_jaws(self):
        if self._optics_is_ready:
            if self._jaws_manually_set:
                jaw_LU = self.jaw_LU
                jaw_RU = self.jaw_RU
                jaw_LD = self.jaw_LD
                jaw_RD = self.jaw_RD
                # TODO: upate gap
                pass
            else:
                # Default to parking
                # TODO: parking is wrt CO, need to correct
                jaw_LU = self.parking
                jaw_RU = -self.parking
                jaw_LD = self.parking
                jaw_RD = -self.parking

                if self.gap_L is not None or self.gap_R is not None:
                    # Get the beam size to be used, depending on align_to
                    if self.align_to == 'maximum':
                        # Reset align_to to the location of the maximum
                        self.align_to = ['front', 'center', 'back'][np.argmax(self.beam_size)]
                    if self.align_to == 'angular':
                        sigma_U = self.beam_size[0]
                        sigma_D = self.beam_size[2]
                    elif self.align_to == 'front':
                        sigma_U = self.beam_size[0]
                        sigma_D = self.beam_size[0]
                    elif self.align_to == 'center':
                        sigma_U = self.beam_size[1]
                        sigma_D = self.beam_size[1]
                    elif self.align_to == 'back':
                        sigma_U = self.beam_size[2]
                        sigma_D = self.beam_size[2]

                    # Get the shift due to the tilt to be used, depending on align_to
                    # TODO: is this correct?
                    # TODO: better name for align_to?
                    if align_to == front:
                        scale = 0
                    elif align_to == back:
                        scale = 1
                    else:
                        scale = 0.5
                    ts_LU = -np.tan(self.tilt_L)*self.active_length*scale
                    ts_RU = -np.tan(self.tilt_R)*self.active_length*scale
                    ts_LD = np.tan(self.tilt_L)*self.active_length*(1-scale)
                    ts_RD = np.tan(self.tilt_R)*self.active_length*(1-scale)

                    if self.onesided in ['left', 'both']:
                        jaw_LU = self.gap_L*sigma_U[0]  + self.offset + ts_LU
                        jaw_LD = self.gap_L*sigma_D[0]  + self.offset + ts_LD
                    if self.onesided in ['right', 'both']:
                        jaw_RU = -self.gap_R*sigma_U[1] + self.offset + ts_RU
                        jaw_RD = -self.gap_R*sigma_D[1] + self.offset + ts_RD

            if self.onesided == 'left':
                self.jaw_LU = min(jaw_LU, self.parking)
                self.jaw_LD = min(jaw_LD, self.parking)
                self.jaw_RU = None
                self.jaw_RD = None
            elif self.onesided == 'right':
                self.jaw_LU = None
                self.jaw_LD = None
                self.jaw_RU = max(jaw_RU, -self.parking)
                self.jaw_RD = max(jaw_RD, -self.parking)
            else:
                self.jaw_LU = min(jaw_LU, self.parking)
                self.jaw_LD = min(jaw_LD, self.parking)
                self.jaw_RU = max(jaw_RU, -self.parking)
                self.jaw_RD = max(jaw_RD, -self.parking)


# Helper functions to set/get properties
# --------------------------------------

def _get_LR(obj, prop, neg=False, name_L='_L', name_R='_R'):
    # Is the property reflected along left/right?
    sign = -1 if neg else 1
    # Is it a property or a dict key?
    if isinstance(obj, dict):
        L = obj[prop + name_L]
        R = obj[prop + name_R]
    else:
        L = getattr(obj, prop + name_L)
        R = getattr(obj, prop + name_R)
    # Find out how many values to return
    if L is None and R is None:
        return None
    elif L is None:
        return R
    elif R is None:
        return L
    else:
        return L if L==sign*R else [L,R]

def _set_LR(obj, prop, val, neg=False, name=None, name_L='_L', name_R='_R'):
    # 'name' is only used for error reporting
    if isinstance(obj, dict):
        name = 'dict_property' if name is None else name
    else:
        name = obj.name if name is None else name
    # Is the property reflected along left/right?
    sign = -1 if neg else 1
    # Find out how to set values
    if not hasattr(val, '__iter__'):
        val = [val]
    if isinstance(val, str):
        raise ValueError(f"Error in settings for {name}: "
                       + f"The setting `{prop}` has to be a number!")
    elif len(val) == 2:
    # The value is of the form [val_L,val_R]
        val_L = val[0]
        val_R = val[1]
    elif len(val) == 1:
    # The value is of the form [val]
        val_L = val[0]
        val_R = val[0]
    else:
        raise ValueError(f"Error in settings for {name}: "
                       + f"The setting `{prop}` must have one or two (L, R) values!")
    # Is it a property or a dict key?
    if isinstance(obj, dict):
        obj[prop + name_L] = val_L
        obj[prop + name_R] = val_R
    else:
        _prop_fset(obj, prop + name_L)(obj, val_L)
        _prop_fset(obj, prop + name_R)(obj, val_R)

def _get_LRUD(obj, prop, neg=False, name=None,
              name_LU='_LU', name_RU='_RU', name_LD='_LD', name_RD='_RD'):
    # 'name' is only used for error reporting
    if isinstance(obj, dict):
        name = 'dict_property' if name is None else name
    else:
        name = obj.name if name is None else name
    # Is the property reflected along left/right?
    sign = -1 if neg else 1
    # Is it a property or a dict key?
    if isinstance(obj, dict):
        LU = obj[prop + name_LU]
        RU = obj[prop + name_RU]
        LD = obj[prop + name_LD]
        RD = obj[prop + name_RD]
    else:
        LU = getattr(obj, prop + name_LU)
        RU = getattr(obj, prop + name_RU)
        LD = getattr(obj, prop + name_LD)
        RD = getattr(obj, prop + name_RD)
    # Find out how many values to return
    if LU is None and RU is None \
    and LD is None and RD is None:
        return None
    elif LU is None and RU is None:
        return LD if LD==sign*RD else [LD,RD]
    elif LD is None and RD is None:
        return LU if LU==sign*RU else [LU,RU]
    elif LU is None and RD is None:
        raise ValueError(f"Error in settings for {name}: "
                       + f"The setting `{prop}` has values for LD and RU but not for LU and RD."
                       + f"Cannot mix jaws L/R with U/D! Either set all four, or only L or only R.")
    elif LD is None and RU is None:
        raise ValueError(f"Error in settings for {name}: "
                       + f"The setting `{prop}` has values for LU and RD but not for LD and RU."
                       + f"Cannot mix jaws L/R with U/D! Either set all four, or only L or only R.")
    else:
        if LU == sign*RU == LD == sign*RD:
            return LU
        elif LU == LD and RU == RD:
            return [LU, RU]
        else:
            return [[LU, RU], [LD, RD]]

def _set_LRUD(obj, prop, val, neg=False, name=None,
              name_LU='_LU', name_RU='_RU', name_LD='_LD', name_RD='_RD'):
    # 'name' is only used for error reporting
    if isinstance(obj, dict):
        name = 'dict_property' if name is None else name
    else:
        name = obj.name if name is None else name
    # Is the property reflected along left/right?
    sign = -1 if neg else 1
    # Find out how to set values
    if not hasattr(val, '__iter__'):
        val = [val]
    if isinstance(val, str):
        raise ValueError(f"Error in settings for {name}: "
                       + f"The setting `{prop}` has to be a number!")
    elif len(val) == 4:
    # The value is of the form [val_LU,val_RU,val_LD,val_RD]
        val_LU = val[0]
        val_RU = val[1]
        val_LD = val[2]
        val_RD = val[3]
    elif len(val) == 2:
        if not hasattr(val[0], '__iter__') and not hasattr(val[1], '__iter__'):
        # The value is of the form [val_L,val_R]
            val_LU = val[0]
            val_RU = val[1]
            val_LD = val[0]
            val_RD = val[1]
        elif not hasattr(val[0], '__iter__') or not hasattr(val[1], '__iter__'):
            raise ValueError(f"Error in settings for {name}: "
                           + f"The setting `{prop}` has to be given as val, [val], "
                           + f"[val_L, val_R], [[val_LU, val_RU], [val_LD, val_RD]], "
                           + f"or [val_LU, val_RU, val_LD, val_RD]!")
        else:
        # The value is of the form [[val_LU, val_RU], [val_LD, val_RD]]
            if isinstance(val[0], str) or isinstance(val[1], str):
                raise ValueError(f"Error in settings for {name}: "
                               + f"The setting `{prop}` has to be a number or a "
                               + f"list of numbers!")
            val_LU = val[0][0]
            val_RU = val[0][1]
            val_LD = val[1][0]
            val_RD = val[1][1]
    elif len(val) == 1:
    # The value is of the form [val]
        val_LU = val[0]
        val_RU = sign*val[0]
        val_LD = val[0]
        val_RD = sign*val[0]
    else:
        raise ValueError(f"Error in settings for {name}: "
                       + f"The setting `{prop}` has to be given as val, [val], "
                       + f"[val_L, val_R], [[val_LU, val_RU], [val_LD, val_RD]], "
                       + f"or [val_LU, val_RU, val_LD, val_RD]!")
    # Is it a property or a dict key?
    if isinstance(obj, dict):
        obj[prop + name_LU] = val_LU
        obj[prop + name_RU] = val_RU
        obj[prop + name_LD] = val_LD
        obj[prop + name_RD] = val_RD
    else:
        _prop_fset(obj, prop + name_LU)(obj, val_LU)
        _prop_fset(obj, prop + name_RU)(obj, val_RU)
        _prop_fset(obj, prop + name_LD)(obj, val_LD)
        _prop_fset(obj, prop + name_RD)(obj, val_RD)

# These getter and setter functions link each property to the corresponding
# entry in the DataFrame
def _prop_fget(obj, attr_name):
    def fget(obj):
        return obj._colldb.loc[obj.name, attr_name]
    return fget

def _prop_fset(obj, attr_name):
    def fset(obj, value):
        obj._colldb.loc[obj.name, attr_name] = value
        # If we want additional effects on the setter funcions, this
        # can be achieved by defining an fset_prop method
        if hasattr(obj.__class__, 'fset_' + attr_name):
            getattr(obj, 'fset_' + attr_name)(value)
        # TODO: self._compute_jaws()
    return fset






class CollimatorDatabase:
    def __init__(self, *, emit, sixtrack_file=None):
        self._optics = pd.DataFrame(columns=['x', 'px', 'y', 'py', 'betx', 'bety', 'alfx', 'alfy', 'dx', 'dy'])
        if sixtrack_file is not None:
            self.load_SixTrack(sixtrack_file)
        else:
            self._colldb = None
        self.emittance = emit
        self._beta_gamma_rel = None

    def __getitem__(self, name):
        return CollimatorSettings(name, self._colldb)

    def to_pandas(self):
        return pd.DataFrame({
                's_center':        self.s_center,
                'gap':             self.gap,
                'jaw':             self.jaw,
                'beam_size':       self.beam_size,
                'aligned_to':      self.align_to,
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


    # TODO: - VALIDATION OF TYPES (e.g. material, stage, align, ..)
    #       - IMPLEMENTATION OF TILT
    #       - CRYSTAL PROPERTIES: only valid if crystal == True (add mask to _set_property?)
    #         make second dataframe for crystals
    #       - show as __repr__

    # The CollDB class has the following fields (those marked
    # with an * are set automatically and cannot be overwritten):
    #   - name
    #   - gap
    #   - jaw *
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

    @property
    def angle(self):
#         angles = np.array([self._colldb.angle_L.values,self._colldb.angle_R.values])
#         return pd.Series([ L if L == R else [L,R] for L, R in angles.T ], index=self._colldb.index, dtype=object)
        return self._colldb['angle_L']

    @angle.setter
    def angle(self, angle):
        self._set_property_LR('angle', angle)
        self._compute_jaws()

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
        self._set_property('offset', offset, single_default_allowed=True)
        self._compute_jaws()

    @property
    def tilt(self):
        tilts = np.array([self._colldb.tilt_L.values,self._colldb.tilt_R.values])
        return pd.Series([ L if L == R else [L,R] for L, R in tilts.T ], index=self._colldb.index, dtype=object)

    @tilt.setter
    def tilt(self, tilts):
        self._set_property_LR('tilt', tilts)
        self._compute_jaws()

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
        self._set_property('parking', parking, single_default_allowed=True)
        self._compute_jaws()

    @property
    def is_active(self):
        return self._colldb['is_active']

    @is_active.setter
    def is_active(self, is_active):
        self._set_property('is_active', is_active, single_default_allowed=True)

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
    def active_length(self, length):
        self._set_property('active_length', length)
        self.align_to = {}

    @property
    def inactive_front(self):
        return self._colldb['inactive_front']

    @inactive_front.setter
    def inactive_front(self, length):
        self._set_property('inactive_front', length)

    @property
    def inactive_back(self):
        return self._colldb['inactive_back']

    @inactive_back.setter
    def inactive_back(self, length):
        self._set_property('inactive_back', length)

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
                raise ValueError("The variable 'gaps' has a different length than the number of collimators in the CollDB. "
                                + "Use a dictionary instead.")
            # Some of the gaps are list (e.g. two different values for both gaps): loop over gaps as dict
            if any(hasattr(gap, '__iter__') for gap in gaps):
                gaps = dict(zip(self.name, gaps))
            # All gaps are single values: use pandas-style assignment
            else:
                # mask those that have an active side for the gap under consideration
                # and have a setting less than 900; the others are set to None
                mask_L = np.logical_and(df.onesided.isin(['both','left']), ~(gaps >= 900))
                mask_R = np.logical_and(df.onesided.isin(['both','right']), ~(gaps >= 900))
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
                gap_L = None if (gap_L is not None and gap_L >= 900) else gap_L
                gap_R = None if (gap_R is not None and gap_R >= 900) else gap_R
                df.loc[name, 'gap_L'] = gap_L
                df.loc[name, 'gap_R'] = gap_R

        if not correct_format:
            raise ValueError("Variable 'gaps' needs to be a pandas Series, dict, numpy array, or list!")

        df.gap_L = df.gap_L.astype('object', copy=False)
        df.gap_R = df.gap_R.astype('object', copy=False)
        self._compute_jaws()

    @property
    def jaw(self):
        jaws = list(np.array([
                        self._colldb.jaw_LU.values,
                        self._colldb.jaw_RU.values,
                        self._colldb.jaw_LD.values,
                        self._colldb.jaw_RD.values
                    ]).T)
        # Need special treatment if there are None's
        def flip(jaw):
            return None if jaw is None else -jaw
        for i, jaw in enumerate(jaws):
            # All 4 jaw points are the same
            if jaw[0] == flip(jaw[1]) == jaw[2] == flip(jaw[3]):
                jaws[i] = jaw[0]
            # Upstream and downstream jaws are the same
            # (all cases except angular alignment and/or tilt)
            elif jaw[0] == jaw[2] and jaw[1] == jaw[3]:
                jaws[i] = [ jaw[0], jaw[1] ]
            else:
                jaws[i] = [ [jaw[0],jaw[1]], [jaw[2],jaw[3]] ]
        return pd.Series(jaws, index=self._colldb.index, dtype=object)

    @property
    def onesided(self):
        return self._colldb.onesided

    @onesided.setter
    def onesided(self, sides):
        self._set_property('onesided', sides, single_default_allowed=True)
        self.gap = self.gap

    @property
    def gamma_rel(self):
        return np.sqrt(self._beta_gamma_rel**2+1)

    @gamma_rel.setter
    def gamma_rel(self, gamma_rel):
        self._beta_gamma_rel = np.sqrt(gamma_rel**2-1)
        self._compute_jaws()

    @property
    def emittance(self):
        return [self._emitx, self._emity]

    @emittance.setter
    def emittance(self, emit):
        if hasattr(emit, '__iter__'):
            if isinstance(emit, str):
                raise ValueError(f"The 'emit' setting has to be a number!")
            elif len(emit) == 2:
                self._emitx = emit[0]
                self._emity = emit[1]
            elif len(emit) == 1:
                self._emitx = emit[0]
                self._emity = emit[0]
            else:
                raise ValueError(f"The 'emit' setting must have one or two values (for emitx and emity)!")
        else:
            self._emitx = emit
            self._emity = emit
        self._compute_jaws()

    @property
    def align_to(self):
        return self._colldb.align_to

    @align_to.setter
    def align_to(self, align):
        self._set_property('align_to', align, single_default_allowed=True, limit_to=['front', 'center', 'back', 'angular'])
        if np.any(self.align_to == 'maximum'):
            raise NotImplementedError
        s_front = self.s_center - self.active_length/2
        s_center = self.s_center
        s_back = self.s_center + self.active_length/2
        mask = self.align_to == 'front'
        self._colldb.loc[mask,'s_align_front'] = s_front[mask]
        self._colldb.loc[mask,'s_align_back']  = s_front[mask]
        mask = self.align_to == 'center'
        self._colldb.loc[mask,'s_align_front'] = s_center[mask]
        self._colldb.loc[mask,'s_align_back']  = s_center[mask]
        mask = self.align_to == 'back'
        self._colldb.loc[mask,'s_align_front'] = s_back[mask]
        self._colldb.loc[mask,'s_align_back']  = s_back[mask]
        mask = self.align_to == 'angular'
        self._colldb.loc[mask,'s_align_front'] = s_front[mask]
        self._colldb.loc[mask,'s_align_back']  = s_back[mask]
        self._compute_jaws()

    # TODO: when does this need to be unset?
    @property
    def _optics_is_ready(self):
        pos = set(self._colldb.s_align_front.values) | set(self._colldb.s_align_back.values)
        return np.all([s in self._optics.index for s in pos]) and self._beta_gamma_rel is not None

    @property
    def betx(self):
        vals = np.array([
            [ self._optics.loc[s,'betx'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'betx'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)

    @property
    def bety(self):
        vals = np.array([
            [ self._optics.loc[s,'bety'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'bety'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)

    @property
    def alfx(self):
        vals = np.array([
            [ self._optics.loc[s,'alfx'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'alfx'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)   

    @property
    def alfy(self):
        vals = np.array([
            [ self._optics.loc[s,'alfy'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'alfy'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)
    
    @property
    def dx(self):
        vals = np.array([
            [ self._optics.loc[s,'dx'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'dx'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)

    @property
    def dy(self):
        vals = np.array([
            [ self._optics.loc[s,'dy'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'dy'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)
    
    @property
    def x(self):
        vals = np.array([
            [ self._optics.loc[s,'x'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'x'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)

    @property
    def px(self):
        vals = np.array([
            [ self._optics.loc[s,'px'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'px'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)

    @property
    def y(self):
        vals = np.array([
            [ self._optics.loc[s,'y'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'y'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)

    @property
    def py(self):
        vals = np.array([
            [ self._optics.loc[s,'py'] if s in self._optics.index else None for s in self._colldb.s_align_front.values ],
            [ self._optics.loc[s,'py'] if s in self._optics.index else None for s in self._colldb.s_align_back.values ]
        ])
        return pd.Series([ F if F == B else [F,B] for F,B in vals.T ], index=self._colldb.index, dtype=object)

    @property
    def beam_size(self):
        if self._optics_is_ready:
            beam_size = np.array([self._beam_size_front,self._beam_size_back])
            return pd.Series([ F if F == B else [F,B] for F, B in beam_size.T ], index=self._colldb.index, dtype=object)
        else:
            return None

    @property
    def _beam_size_front(self):
        df = self._colldb
        opt = self._optics
        betx = opt.loc[df.s_align_front,'betx'].astype(float)
        bety = opt.loc[df.s_align_front,'bety'].astype(float)
        sigmax = np.sqrt(betx*self._emitx/self._beta_gamma_rel)
        sigmay = np.sqrt(bety*self._emity/self._beta_gamma_rel)
        result = np.sqrt(
                    (sigmax*np.cos(np.float_(df.angle.values)*np.pi/180))**2
                    + (sigmay*np.sin(np.float_(df.angle.values)*np.pi/180))**2
                )
        result.index = self._colldb.index
        return result

    @property
    def _beam_size_back(self):
        df = self._colldb
        opt = self._optics
        betx = opt.loc[df.s_align_back,'betx'].astype(float)
        bety = opt.loc[df.s_align_back,'bety'].astype(float)
        sigmax = np.sqrt(betx*self._emitx/self._beta_gamma_rel)
        sigmay = np.sqrt(bety*self._emity/self._beta_gamma_rel)
        result = np.sqrt(
                    (sigmax*np.cos(np.float_(df.angle.values)*np.pi/180))**2
                    + (sigmay*np.sin(np.float_(df.angle.values)*np.pi/180))**2
                )
        result.index = self._colldb.index
        return result

    # parking is defined with respect to closed orbit
    # TODO: tilt
    # 'upstr'  =>  'front'  en   'downstr'  =>  'back'
    def _compute_jaws(self):
        if self._optics_is_ready:
            df = self._colldb
            beam_size_front = self._beam_size_front
            beam_size_back  = self._beam_size_back
            jaw_LU = df['gap_L']*beam_size_front + self.offset
            jaw_RU = df['gap_R']*beam_size_front - self.offset
            jaw_LD = df['gap_L']*beam_size_back  + self.offset
            jaw_RD = df['gap_R']*beam_size_back  - self.offset
            df['jaw_LU'] = df['parking'] if df['gap_L'] is None else np.minimum(jaw_LU,df['parking'])
            df['jaw_RU'] = -df['parking'] if df['gap_R'] is None else -np.minimum(jaw_RU,df['parking'])
            df['jaw_LD'] = df['parking'] if df['gap_L'] is None else np.minimum(jaw_LD,df['parking'])
            df['jaw_RD'] = -df['parking'] if df['gap_R'] is None else -np.minimum(jaw_RD,df['parking'])



    # ---------------------------------------
    # ------ Property setter functions ------
    # ---------------------------------------

    def _set_property(self, prop, vals, single_default_allowed=False, limit_to=[]):
        df = self._colldb
        if not isinstance(limit_to, (list, tuple, set)):
            limit_to = [limit_to]
        if isinstance(vals, dict):
            for name, val in vals.items():
                if name not in self.name:
                    raise ValueError(f"Collimator {name} not found in CollDB!")
                if limit_to!=[] and val not in limit_to:
                    raise ValueError(f"Cannot set {prop} to {val}. Choose from {limit_to}!")
                df.loc[name, prop] = val
        elif isinstance(vals, pd.Series) or isinstance(vals, list) or isinstance(vals, np.ndarray):
            if len(vals) != len(self.name):
                raise ValueError(f"The variable '{prop}' has a different length than the number of collimators in the CollDB. "
                                + "Use a dictionary instead.")
            if limit_to!=[] and np.any([val not in limit_to for val in vals]):
                raise ValueError(f"Cannot set {prop} to {vals}. Choose from {limit_to}!")
            df[prop] = vals
        else:
            if single_default_allowed:
                if limit_to!=[] and vals not in limit_to:
                    raise ValueError(f"Cannot set {prop} to {vals}. Choose from {limit_to}!")
                df[prop] = vals
            else:
                raise ValueError(f"Variable '{prop}' needs to be a pandas Series, dict, numpy array, or list!")


    def _set_property_LR(self, prop, vals):
        df = self._colldb
        correct_format = False
        # The variable vals is a Series or a list
        if isinstance(vals, pd.Series) or isinstance(vals, list) or isinstance(vals, np.ndarray):
            correct_format = True
            if len(vals) != len(self.name):
                raise ValueError(f"The variable '{prop}' has a different length than the number of collimators in the CollDB. "
                                + "Use a dictionary instead.")
            # Some of the vals are list (e.g. two different values for both gaps): loop over vals as dict
            if any(hasattr(val, '__iter__') for val in vals):
                vals = dict(zip(self.name, vals))
            # All gaps are single values: use pandas-style assignment
            else:
                df[prop + "_L"] = vals
                df[prop + "_R"] = vals
        
        # The variable vals is a dictionary
        if isinstance(vals, dict):
            correct_format = True
            for name, val in vals.items():
                if name not in self.name:
                    raise ValueError(f"Collimator {name} not found in CollDB!")
                if hasattr(val, '__iter__'):
                    if isinstance(val, str):
                        raise ValueError(f"The '{prop}' setting has to be a number!")
                    elif len(val) == 2:
                        val_L = val[0]
                        val_R = val[1]
                    elif len(val) == 1:
                        val_L = val[0]
                        val_R = val[0]
                    else:
                        raise ValueError(f"The '{prop}' setting must have one or two values (for the left and the right jaw)!")
                else:
                    val_L = val
                    val_R = val
                df.loc[name, prop + "_L"] = val_L
                df.loc[name, prop + "_R"] = val_R

        if not correct_format:
            raise ValueError("Variable '{prop}' needs to be a pandas Series, dict, numpy array, or list!")

        df[prop + "_L"] = df[prop + "_L"].astype('object', copy=False)
        df[prop + "_R"] = df[prop + "_R"].astype('object', copy=False)
                
            
    def _initialise_None(self):
        fields = {'s_center':None, 'align_to': None, 's_align_front': None, 's_align_back': None }
        fields.update({'gap_L': None, 'gap_R': None, 'angle_L': 0, 'angle_R': 0, 'offset': 0, 'tilt_L': 0, 'tilt_R': 0})
        fields.update({'jaw_LU': None, 'jaw_RU': None, 'jaw_LD': None, 'jaw_RD': None, 'parking': None, 'family': None})
        fields.update({'onesided': 'both', 'material': None, 'stage': None, 'collimator_type': None, 'is_active': True})
        fields.update({'active_length': 0, 'inactive_front': 0, 'inactive_back': 0, 'sigmax': None, 'sigmay': None})
        fields.update({'crystal': None, 'bend': None, 'xdim': 0, 'ydim': 0, 'miscut': 0, 'thick': 0})
        fields.update({'ref_xU': 0, 'ref_yU': 0, 'ref_xD': 0, 'ref_yD': 0})
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
                if len(sline) > 0 and len(sline) < 6:
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

        names = ['name', 'jaw', 'material', 'length', 'angle', 'offset']

        df = pd.read_csv(io.StringIO(coll_data_string), delim_whitespace=True,
                        index_col=False, names=names)

        df = df[['name', 'jaw', 'length', 'angle', 'material', 'offset']]
        df.insert(5,'stage', df['jaw'].apply(lambda s: family_types.get(s, 'UNKNOWN')))   

        sides = df['name'].apply(lambda s: onesided.get(s, 0))
        gaps = df['jaw'].apply(lambda s: float(family_settings.get(s, s)))

        df['name'] = df['name'].str.lower() # Make the names lowercase for easy processing
        df['angle_L'] = df['angle']
        df['angle_R'] = df['angle']
        df.rename(columns={'length':'active_length', 'offset':'offset'}, inplace=True)
        df['parking'] = 0.025
        df.loc[df.name.str[:3] == 'tct', 'parking'] = 0.04

        df = df.set_index('name')
        self._colldb = df.drop('jaw', axis=1)

        self._initialise_None()
        self.gap = gaps.values
        self._colldb.onesided = sides.values
        self._colldb.onesided = [ 'both' if s==0 else s for s in self._colldb.onesided ]
        self._colldb.onesided = [ 'left' if s==1 else s for s in self._colldb.onesided ]
        self._colldb.onesided = [ 'right' if s==2 else s for s in self._colldb.onesided ]
        
        self.gap = self.gap

        # Check if collimator active
        # Check if gap is list (assymetric jaws)
