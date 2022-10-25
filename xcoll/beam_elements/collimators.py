import numpy as np
# from abc import ABC, ABCMeta, abstractmethod

import xobjects as xo
import xtrack as xt

from ..general import _pkg_root
from ..tables import CollimatorImpacts

# class MetaCollimator(xt.base_element.MetaBeamElement, ABCMeta):
#     pass

class BaseCollimator(xt.BeamElement):#, metaclass=MetaCollimator):
    _xofields = {
        'inactive_front': xo.Float64,
        'active_length': xo.Float64,
        'inactive_back': xo.Float64,
        'jaw_F_L': xo.Float64,
        'jaw_F_R': xo.Float64,
        'jaw_B_L': xo.Float64,
        'jaw_B_R': xo.Float64,
        'jaw_U': xo.Float64,
        'jaw_D': xo.Float64,
        'dx': xo.Float64,
        'dy': xo.Float64,
        'cos_z': xo.Float64,
        'sin_z': xo.Float64,
        '_active': xo.Int8,
        '_record_impacts': xo.Int8,
        '_impacts': xo.Ref(CollimatorImpacts._XoStruct)
    }

    isthick = True
    behaves_like_drift = True
    # TODO: how to pass _impacts to from_dict()... ?
    _skip_in_to_dict = ['_impacts', '_active', '_record_impacts']
    _store_in_to_dict = ['angle', 'is_active']

    def __init__(self, **kwargs):
        # TODO: quick hack to avoid instantiation; did not manage to get it to work correclty with ABC
        if self.__class__.__name__ == 'BaseCollimator':
            raise Exception("Abstract class 'BaseCollimator' cannot be instantiated!")
        kwargs.setdefault('jaw_F_L', 1)
        kwargs.setdefault('jaw_F_R', -1)
        kwargs.setdefault('jaw_B_L', 1)
        kwargs.setdefault('jaw_B_R', -1)
        kwargs.setdefault('jaw_U', 1)
        kwargs.setdefault('jaw_D', -1)
        kwargs.setdefault('inactive_front', 0)
        kwargs.setdefault('inactive_back', 0)
        kwargs.setdefault('dx', 0)
        kwargs.setdefault('dy', 0)
        angle = kwargs.pop('angle', 0)
        is_active = kwargs.pop('is_active', True)
        impacts = kwargs.pop('impacts', None)

        anglerad = angle / 180. * np.pi
        kwargs['cos_z'] = np.cos(anglerad)
        kwargs['sin_z'] = np.sin(anglerad)
        is_active = 1 if is_active == True else is_active
        is_active = 0 if is_active == False else is_active
        kwargs['_active'] = is_active
        if impacts is None:
            kwargs['_record_impacts'] = 0
        else:
            kwargs['_record_impacts'] = 1
        kwargs['_impacts'] = impacts
        super().__init__(**kwargs)


    @property
    def angle(self):
        return np.arctan2(self.sin_z, self.cos_z) * (180.0 / np.pi)

    @angle.setter
    def angle(self, angle):
        anglerad = angle / 180. * np.pi
        self.cos_z = np.cos(anglerad)
        self.sin_z = np.sin(anglerad)

    @property
    def is_active(self):
        return True if self._active == 1 else False

    @is_active.setter
    def is_active(self, is_active):
        is_active = 1 if is_active == True else is_active
        is_active = 0 if is_active == False else is_active
        self._active = is_active
        if is_active <= 0:
            self.jaw_F_L = 1
            self.jaw_F_R = -1
            self.jaw_B_L = 1
            self.jaw_B_R = -1 

    @property
    def length(self):
        return (self.inactive_front + self.active_length + self.inactive_back)

    @property
    def impacts(self):
        return self._impacts

    @impacts.setter
    def impacts(self, impacts):
        if impacts is None:
            self._record_impacts = 0
        elif isinstance(impacts, CollimatorImpacts):
            self._record_impacts = 1
        else:
            raise ValueError("The variable 'impacts' needs to be a CollimatorImpacts object!")
        self._impacts = impacts


class BlackAbsorber(BaseCollimator):
    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements/collimators_src/absorber.h')
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)





