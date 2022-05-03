import numpy as np

import xobjects as xo
import xtrack as xt

from ..general import _pkg_root
from ..collimator_impacts import CollimatorImpactsData

class BlackAbsorber(xt.BeamElement):
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
        'impacts': xo.Ref(CollimatorImpactsData)
    }

    isthick = True
    behaves_like_drift = True

    def __init__(self, angle=0, is_active=True, impacts=None,
                 jaw_F_L=1, jaw_F_R=-1, jaw_B_L=1, jaw_B_R=-1, jaw_U=1, jaw_D=-1,
                 **kwargs):
        super().__init__(**kwargs)
        self.angle = angle
        self.jaw_F_L = jaw_F_L
        self.jaw_F_R = jaw_F_R
        self.jaw_B_L = jaw_B_L
        self.jaw_B_R = jaw_B_R
        self.jaw_U   = jaw_U
        self.jaw_D   = jaw_D
        self.is_active = is_active
        if impacts is None:
            self._record_impacts = 0
        else:
            self._record_impacts = 1
        self.impacts = impacts
        

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

BlackAbsorber.XoStruct.extra_sources = [
        _pkg_root.joinpath('beam_elements/collimators_src/absorber.h')]