import numpy as np

import xobjects as xo
import xtrack as xt

from ..general import _pkg_root

class Collimator(xt.BeamElement):
    _xofields = {
        'inactive_length_at_start': xo.Float64,
        'active_length': xo.Float64,
        'inactive_length_at_end': xo.Float64,
        'n_slices': xo.Int64,
        'a_min': xo.Float64,
        'a_max': xo.Float64,
        'b_min': xo.Float64,
        'b_max': xo.Float64,
        'dx': xo.Float64,
        'dy': xo.Float64,
        'cos_z': xo.Float64,
        'sin_z': xo.Float64,
    }

    isthick = True
    behaves_like_drift = True

    def __init__(self, angle=0, **kwargs):
        anglerad = angle / 180 * np.pi
        kwargs['cos_z']=np.cos(anglerad)
        kwargs['sin_z']=np.sin(anglerad)
        super().__init__(**kwargs)

    @property
    def angle(self):
        return np.arctan2(self.sin_z, self.cos_z) * (180.0 / np.pi)

    @property
    def length(self):
        return (self.inactive_length_at_start + self.active_length
                + self.inactive_length_at_start)

Collimator.XoStruct.extra_sources = [
        _pkg_root.joinpath('beam_elements/collimator_src/collimator.h')]