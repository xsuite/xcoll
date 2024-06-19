# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xpart as xp
import xtrack as xt

from .base import BaseCollimator, BaseCrystal, InvalidXcoll
from ..scattering_routines.k2 import K2Engine, track
from ..general import _pkg_root


class _K2Collimator(BaseCollimator):
    _xofields = {**BaseCollimator._xofields,
        '_material':        xo.String,
        '_tracking':        xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [ *BaseCollimator._skip_in_to_dict, '_material' ]
    _store_in_to_dict      = [ *BaseCollimator._store_in_to_dict, 'material' ]
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, K2Engine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','everest_collimator.h')
    ]

    _kernels = {
        'EverestCollimator_set_material': xo.Kernel(
                c_name='EverestCollimator_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = Material()
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)

    @property
    def track(self, part):
        if self._tracking:
            track(self, part)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                                 _buffer=_buffer, _offset=_offset)



# class EverestCrystal(BaseCrystal):
#     _xofields = {**BaseCrystal._xofields,
#         'miscut':             xo.Float64,
#         '_orient':            xo.Int8,
#         '_critical_angle':    xo.Float64,
#         '_material':          CrystalMaterial,
#         'rutherford_rng':     xt.RandomRutherford,
#         '_tracking':          xo.Int8
#     }

#     isthick = True
#     needs_rng = True
#     allow_track = True
#     behaves_like_drift = True
#     skip_in_loss_location_refinement = True

#     _skip_in_to_dict       = [*BaseCrystal._skip_in_to_dict, '_orient', '_material']
#     _store_in_to_dict      = [*BaseCrystal._store_in_to_dict, 'lattice', 'material']
#     _internal_record_class = BaseCrystal._internal_record_class

#     _depends_on = [BaseCrystal, EverestEngine]

#     _extra_c_sources = [
#         _pkg_root.joinpath('beam_elements','elements_src','everest_crystal.h')
#     ]

#     _kernels = {
#         'EverestCrystal_set_material': xo.Kernel(
#                 c_name='EverestCrystal_set_material',
#                 args=[xo.Arg(xo.ThisClass, name='el')]
#             )
#         }


#     def __init__(self, **kwargs):
#         to_assign = {}
#         if '_xobject' not in kwargs:
#             to_assign['material'] = kwargs.pop('material', None)
#             kwargs['_material'] = CrystalMaterial()
#             to_assign['lattice'] = kwargs.pop('lattice', 'strip')
#             kwargs.setdefault('miscut', 0)
#             kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
#             kwargs.setdefault('_tracking', True)
#         super().__init__(**kwargs)
#         for key, val in to_assign.items():
#             setattr(self, key, val)


#     @property
#     def material(self):
#         return self._material

#     @material.setter
#     def material(self, material):
#         if material is None:
#             material = CrystalMaterial()
#         if isinstance(material, dict):
#             material = CrystalMaterial.from_dict(material)
#         if not isinstance(material, CrystalMaterial):
#             raise ValueError("Invalid material!")
#         if not xt.line._dicts_equal(self.material.to_dict(), material.to_dict()):
#             self._material = material
#             self.EverestCrystal_set_material(el=self)

#     @property
#     def critical_angle(self):
#         return self._critical_angle if abs(self._critical_angle) > 1.e-10 else None

#     @property
#     def lattice(self):
#         if self._orient == 1:
#             return 'strip'
#         elif self._orient == 2:
#             return 'quasi-mosaic'
#         else:
#             raise ValueError(f"Illegal value {self._orient} for '_orient'!")

#     @lattice.setter
#     def lattice(self, lattice):
#         if lattice == 'strip' or lattice == '110' or lattice == 110:
#             self._orient = 1
#         elif lattice == 'quasi-mosaic' or lattice == '111' or lattice == 111:
#             self._orient = 2
#         else:
#             raise ValueError(f"Illegal value {lattice} for 'lattice'! "
#                             + "Only use 'strip' (110) or 'quasi-mosaic' (111).")


#     def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
#         return InvalidXcoll(length=-self.length, _context=_context,
#                                  _buffer=_buffer, _offset=_offset)


