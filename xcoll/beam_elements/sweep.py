# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt
from xtrack import BeamElement, Cavity
from xtrack.beam_elements.elements import _HasModelRF, _HasIntegrator

from ..general import _pkg_root


class SweepCavity(_HasModelRF, _HasIntegrator, BeamElement):
    _xofields = {
        'length':      xo.Float64,
        'cavity':      Cavity._XoStruct,
        'df_per_turn': xo.Float64,     # Fixed df per turn
        'df':          xo.Float64[:],  # df at each turn, for turn-dependent frequency control
        '_last_turn':  xo.Int64,
    }

    isthick = True
    has_backtrack = True
    allow_loss_refinement = True

    _noexpr_fields   = Cavity._noexpr_fields
    _skip_in_to_dict = Cavity._skip_in_to_dict

    _depends_on = [Cavity]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','sweep_cavity.h')
    ]

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            if 'cavity' not in kwargs:
                raise ValueError("Need to provide an existing cavity.")
            if not isinstance(kwargs['cavity'], Cavity):
                raise ValueError("`cavity` must be an instance of xtrack.Cavity!")
            if 'length' in kwargs and kwargs['length'] != kwargs['cavity'].length:
                raise ValueError("Length of SweepCavity must be the same as that of the provided cavity.")
            kwargs['length'] = kwargs['cavity'].length
            if 'df_per_turn' in kwargs:
                if 'df' in kwargs:
                    raise ValueError("Provide either `df` or `df_per_turn`, not both.")
                kwargs['_last_turn'] = - 1
                kwargs['df'] = np.array([])
            elif 'df' in kwargs:
                if not hasattr(kwargs['df'], '__iter__') and not isinstance(kwargs['df'], str):
                    kwargs['df'] = [kwargs['df']]
                kwargs['df'] = np.asarray(kwargs['df'], dtype=np.float64)
                kwargs['_last_turn'] = len(kwargs['df']) - 1
            else:
                raise ValueError("Need to provide `df` or `df_per_turn`.")
        super().__init__(**kwargs)

    def __getattr__(self, name):
        if hasattr(self.cavity, name):
            return getattr(self.cavity, name)
        else:
            return super().__getattr__(name)