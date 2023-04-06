# copyright ############################### #
# This file is part of the Xcoll Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import xobjects as xo
import xpart as xp
import xtrack as xt

from .base_collimator import BaseCollimator, InvalidCollimator
from ..scattering_routines.everest import Material, CrystalMaterial, EverestEngine
from ..general import _pkg_root



# TODO: 
#    - remove dx, dy, offset, tilt, as this should only be in colldb (and here only the jaw positions)
#    - We want these elements to behave as if 'iscollective = True' when doing twiss etc (because they would ruin the CO),
#      but as if 'iscollective = False' for normal tracking as it is natively in C...
#      Currently this is achieved with the hack '_tracking' which defaults to False after installation in the line, and is
#      only activated around the track command. Furthermore, because of 'iscollective = False' we need to specify get_backtrack_element
#      We want it nicer..

class EverestCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'offset':           xo.Float64,
        'onesided':         xo.Int8,
        'tilt':             xo.Float64[:],  # TODO: how to limit this to length 2
        'material':         Material,
        'rutherford_rng':   xt.RandomRutherford,
        '_tracking':        xo.Int8
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

#     _skip_in_to_dict       = [ *BaseCollimator._skip_in_to_dict, '_material' ]
#     _store_in_to_dict      = [ *BaseCollimator._store_in_to_dict, 'material' ]
    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','everest_collimator.h')
    ]


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('offset', 0)
            kwargs.setdefault('onesided', False)
            kwargs.setdefault('tilt', [0,0])
            tilt = kwargs['tilt']
            if hasattr(tilt, '__iter__'):
                if isinstance(tilt, str):
                    raise ValueError("Variable tilt has to be a number or array of numbers!")
                elif len(tilt) == 1:
                    tilt = [tilt[0], tilt[0]]
                elif len(tilt) > 2:
                    raise ValueError("Variable tilt cannot have more than two elements (tilt_L and tilt_R)!")
            else:
                tilt = [tilt, tilt]
            kwargs['tilt'] = tilt
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
#         if '_xobject' not in kwargs:
#             self.random_generator.set_rutherford_by_xcoll_material(self.material)

#     @property
#     def material(self):
#         return self._material

#     @material.setter
#     def material(self, material):
#         self._material = material
# #         self.random_generator.set_rutherford_by_xcoll_material(material)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        # TODO: this should be an InvalidCollimator
        return xt.Drift(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)



class EverestCrystal(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'align_angle':    xo.Float64,  #  = - sqrt(eps/beta)*alpha*nsigma
        'bend':           xo.Float64,
        'xdim':           xo.Float64,
        'ydim':           xo.Float64,
        'thick':          xo.Float64,
        'crytilt':        xo.Float64,
        'miscut':         xo.Float64,
        'orient':         xo.Float64,
        'offset':         xo.Float64,
        'onesided':       xo.Int8,
        'tilt':           xo.Float64[:],  # TODO: how to limit this to length 2
        'material':       CrystalMaterial,
        'rutherford_rng': xt.RandomRutherford,
        '_tracking':      xo.Int8
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

#     _skip_in_to_dict       = [ *BaseCollimator._skip_in_to_dict, '_material' ]
#     _store_in_to_dict      = [ *BaseCollimator._store_in_to_dict, 'material' ]
    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','everest_crystal.h')
    ]


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('offset', 0)
            kwargs.setdefault('onesided', False)
            kwargs.setdefault('tilt', [0,0])
            tilt = kwargs['tilt']
            if hasattr(tilt, '__iter__'):
                if isinstance(tilt, str):
                    raise ValueError("Variable tilt has to be a number or array of numbers!")
                elif len(tilt) == 1:
                    tilt = [tilt[0], tilt[0]]
                elif len(tilt) > 2:
                    raise ValueError("Variable tilt cannot have more than two elements (tilt_L and tilt_R)!")
            else:
                tilt = [tilt, tilt]
            kwargs['tilt'] = tilt
            kwargs.setdefault('bend', 0)
            kwargs.setdefault('xdim', 0)
            kwargs.setdefault('ydim', 0)
            kwargs.setdefault('thick', 0)
            kwargs.setdefault('crytilt', 0)
            kwargs.setdefault('miscut', 0)
            kwargs.setdefault('orient', 0)
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
#         if '_xobject' not in kwargs:
#             self.random_generator.set_rutherford_by_xcoll_material(self.material)

#     @property
#     def material(self):
#         return self._material

#     @material.setter
#     def material(self, material):
#         self._material = material
#         self.random_generator.set_rutherford_by_xcoll_material(material)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        # TODO: this should be an InvalidCollimator
        return xt.Drift(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)


