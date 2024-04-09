import io
import json
import numpy as np
import pandas as pd


# in colldb:
#      (gap [sigma] + offset [m] + tilt [deg]) or opening [m]
#      angle
#      length (all three)
#      side
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
#      side
#      material
#      is_active

# in neither (as they are sequence/session-dependent):
#      align_to
#      s_center
#      collimator_type

# if physical_opening is used, other variables like tilt, opening, offset, .. are ignored
_element_properties  = {
                        'length':           0,
                        'jaw':              None,
                        'angle':            0,
                        'side':             'both',
                        'material':         None,
                        'active':           True,
                        'crystal':          False
                       }
_crystal_properties  = {
                        'bend':             None,
                        'xdim':             0,
                        'ydim':             0,
                        'miscut':           0,
                        'thick':            0,
#         'align_angle':    xo.Float64,  #  = - sqrt(eps/beta)*alpha*nsigma
#         'crytilt':        xo.Float64,
#         'orient':         xo.Float64,
                       }
_sequence_properties = {
                        's_center':         None,
                        'align_to':         None,
                        'collimator_type':  None,
                       }
_colldb_properties   = {
                        'gap':              None,  # [L, R]
                        'offset':           0,     # single value: shift of gap in mm
                        'extra_mm':         0,     # single value: widening of gap
                        'opening_mm':       None,  # [L, R]
                        'tilt':             0,
                        'parking':          1,
                        'stage':            None,
                        'family':           None,
                        'overwritten_keys': []
                       }
_optics_vals = ['x', 'px', 'y', 'py', 'betx', 'bety', 'alfx', 'alfy', 'dx', 'dy', 'sigma_x', 'sigma_y']


# This creates a view on the settings of one collimator, kept in sync with the main database
class CollimatorSettings:

    def __init__(self, name, colldb=None, optics=None, element=None):
        if colldb is None:
            colldb = pd.DataFrame({**_element_properties, **_crystal_properties}, index=[name])
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
                    ts_LU = -np.tan(self.tilt_L)*self.length*scale
                    ts_RU = -np.tan(self.tilt_R)*self.length*scale
                    ts_LD = np.tan(self.tilt_L)*self.length*(1-scale)
                    ts_RD = np.tan(self.tilt_R)*self.length*(1-scale)

                    if self.side in ['left', 'both']:
                        jaw_LU = self.gap_L*sigma_U[0]  + self.offset + ts_LU
                        jaw_LD = self.gap_L*sigma_D[0]  + self.offset + ts_LD
                    if self.side in ['right', 'both']:
                        jaw_RU = -self.gap_R*sigma_U[1] + self.offset + ts_RU
                        jaw_RD = -self.gap_R*sigma_D[1] + self.offset + ts_RD

            if self.side == 'left':
                self.jaw_LU = min(jaw_LU, self.parking)
                self.jaw_LD = min(jaw_LD, self.parking)
                self.jaw_RU = None
                self.jaw_RD = None
            elif self.side == 'right':
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
    if name is None:
        if isinstance(obj, dict):
            name = 'dict_property'
        else:
            name = obj.name if hasattr(obj, 'name') else obj.__class__.__name__
    # Is the property reflected along left/right?
    sign = -1 if neg else 1
    # Find out how to set values
    if not hasattr(val, '__iter__'):
        val = [val]
    error_dimension = False
    if isinstance(val, str):
        raise ValueError(f"Error in settings for {name}: "
                       + f"The setting `{prop}` has to be a number!")
    elif len(val) == 2:
    # The value is of the form [val_L,val_R]
        if hasattr(val[0], '__iter__') or hasattr(val[1], '__iter__'):
            error_dimension = True
        val_L = val[0]
        val_R = val[1]
    elif len(val) == 1:
    # The value is of the form [val]
        if hasattr(val[0], '__iter__'):
            error_dimension = True
        val_L = val[0]
        val_R = sign*val[0]
    else:
        error_dimension = True
    if error_dimension:
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
    if name is None:
        if isinstance(obj, dict):
            name = 'dict_property'
        else:
            name = obj.name if hasattr(obj, 'name') else obj.__class__.__name__
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
    if name is None:
        if isinstance(obj, dict):
            name = 'dict_property'
        else:
            name = obj.name if hasattr(obj, 'name') else obj.__class__.__name__
    # Is the property reflected along left/right?
    sign = -1 if neg else 1

    # Find out how to set values
    if not hasattr(val, '__iter__'):
        val = [val]
    error_string    = False
    error_dimension = False
    if isinstance(val, str):
        error_string = True
    elif len(val) == 4:
    # The value is of the form [val_LU,val_RU,val_LD,val_RD]
        if hasattr(val[0], '__iter__') or hasattr(val[1], '__iter__') \
        or hasattr(val[2], '__iter__') or hasattr(val[3], '__iter__'):
            error_dimension = True
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
            error_dimension = True
        else:
        # The value is of the form [[val_LU, val_RU], [val_LD, val_RD]]
            if isinstance(val[0], str) or isinstance(val[1], str):
                error_string = True
            if hasattr(val[0][0], '__iter__') or hasattr(val[0][1], '__iter__') \
            or hasattr(val[1][0], '__iter__') or hasattr(val[1][1], '__iter__'):
                error_dimension = True
            val_LU = val[0][0]
            val_RU = val[0][1]
            val_LD = val[1][0]
            val_RD = val[1][1]
    elif len(val) == 1:
    # The value is of the form [val]
        if hasattr(val[0], '__iter__'):
            error_dimension = True
        val_LU = val[0]
        val_RU = sign*val[0]
        val_LD = val[0]
        val_RD = sign*val[0]
    else:
        error_dimension = True
    if error_dimension:
        raise ValueError(f"Error in settings for {name}: "
                       + f"The setting `{prop}` has to be given as val, [val], "
                       + f"[val_L, val_R], [[val_LU, val_RU], [val_LD, val_RD]], "
                       + f"or [val_LU, val_RU, val_LD, val_RD]!")
    if error_string:
        raise ValueError(f"Error in settings for {name}: "
                               + f"The setting `{prop}` has to be a number or a "
                               + f"list of numbers!")

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
        if hasattr(obj, '_colldb'):
            obj._colldb.loc[obj.name, attr_name] = value
        if hasattr(obj, attr_name):
            setattr(obj, attr_name, value)
        # If we want additional effects on the setter funcions, this
        # can be achieved by defining an fset_prop method
        if hasattr(obj.__class__, 'fset_' + attr_name):
            getattr(obj, 'fset_' + attr_name)(value)
        # TODO: self._compute_jaws()
    return fset

