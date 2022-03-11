import xobjects as xo
import xtrack as xt

from ..general import _pkg_root

class Collimator(xt.BeamElement):
    _xofields = {
        'inactive_length_at_start': xo.Float64,
        'active_length': xo.Float64,
        'inactive_length_at_end': xo.Float64,
        'n_slices': xo.Int64,
        'a_max': xo.Float64,
        'a_min': xo.Float64,
        'b_max': xo.Float64,
        'b_min': xo.Float64,
        'cos_z': xo.Float64,
        'sin_z': xo.Float64,
    }

Collimator.XoStruct.extra_sources = [
        _pkg_root.joinpath('beam_elements/collimator_src/collimator.h')]