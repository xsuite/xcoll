# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from .base import BaseBlock, BaseCollimator, BaseCrystal
from ..scattering_routines.everest import Material, CrystalMaterial, EverestEngine
from ..general import _pkg_root


# TODO:
#      We want these elements to behave as if 'iscollective = True' when doing twiss etc (because they would ruin the CO),
#      but as if 'iscollective = False' for normal tracking as it is natively in C...
#      Currently this is achieved with the hack '_tracking' which defaults to False after installation in the line, and is
#      only activated around the track command.


class EverestBlock(BaseBlock):
    _xofields = BaseBlock._xofields | {
        '_material':        Material,
        'rutherford_rng':   xt.RandomRutherford,
        '_tracking':        xo.Int8
    }

    needs_rng         = True
    allow_track       = True
    _depends_on       = [BaseBlock, EverestEngine]
    _noexpr_fields    = BaseBlock._noexpr_fields | {'material'}
    _skip_in_to_dict  = BaseBlock._skip_in_to_dict + ['_material']
    _store_in_to_dict = BaseBlock._store_in_to_dict + ['material']
    _extra_c_sources  = [_pkg_root / 'beam_elements' / 'elements_src' / 'everest_block.h']

    _kernels = {
        'EverestBlock_set_material': xo.Kernel(
                c_name='EverestBlock_set_material',
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
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if material is None:
            material = Material()
        if isinstance(material, dict):
            material = Material.from_dict(material)
        if not isinstance(material, Material):
            raise ValueError("Invalid material!")
        if not xt.line._dicts_equal(self.material.to_dict(), material.to_dict()):
            self._material = material
            self.EverestBlock_set_material(el=self)


class EverestCollimator(BaseCollimator):
    _xofields = EverestBlock._xofields | BaseCollimator._xofields

    needs_rng         = True
    allow_track       = True
    _depends_on       = [BaseCollimator, EverestEngine]
    _noexpr_fields    = BaseCollimator._noexpr_fields | {'material'}
    _skip_in_to_dict  = BaseCollimator._skip_in_to_dict + ['_material']
    _store_in_to_dict = BaseCollimator._store_in_to_dict + ['material']
    _extra_c_sources  = [_pkg_root / 'beam_elements' / 'elements_src' / 'everest_collimator.h']

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
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if material is None:
            material = Material()
        if isinstance(material, dict):
            material = Material.from_dict(material)
        if not isinstance(material, Material):
            raise ValueError("Invalid material!")
        if not xt.line._dicts_equal(self.material.to_dict(), material.to_dict()):
            self._material = material
            self.EverestCollimator_set_material(el=self)


class EverestCrystal(BaseCrystal):
    _xofields = EverestBlock._xofields | BaseCrystal._xofields | {
        '_material':          CrystalMaterial,
        'miscut':             xo.Float64,
        '_orient':            xo.Int8,
        '_critical_angle':    xo.Float64,
        '_critical_radius':   xo.Float64
    }

    needs_rng         = True
    allow_track       = True
    _depends_on       = [BaseCrystal, EverestEngine]
    _noexpr_fields    = BaseCrystal._noexpr_fields | {'material', 'lattice'}
    _skip_in_to_dict  = BaseCrystal._skip_in_to_dict + ['_material', '_orient']
    _store_in_to_dict = BaseCrystal._store_in_to_dict + ['material', 'lattice']
    _extra_c_sources  = [_pkg_root / 'beam_elements' / 'elements_src' / 'everest_crystal.h']

    _kernels = {
        'EverestCrystal_set_material': xo.Kernel(
                c_name='EverestCrystal_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = CrystalMaterial()
            to_assign['lattice'] = kwargs.pop('lattice', 'strip')
            kwargs.setdefault('miscut', 0)
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if material is None:
            material = CrystalMaterial()
        if isinstance(material, dict):
            material = CrystalMaterial.from_dict(material)
        if not isinstance(material, CrystalMaterial):
            raise ValueError("Invalid material!")
        if not xt.line._dicts_equal(self.material.to_dict(), material.to_dict()):
            self._material = material
            self.EverestCrystal_set_material(el=self)

    @property
    def critical_angle(self):
        return self._critical_angle if abs(self._critical_angle) > 1.e-12 else None

    @property
    def critical_radius(self):
        return self._critical_radius if abs(self._critical_radius) > 1.e-10 else None

    @property
    def lattice(self):
        if self._orient == 1:
            return 'strip'
        elif self._orient == 2:
            return 'quasi-mosaic'
        else:
            raise ValueError(f"Illegal value {self._orient} for '_orient'!")

    @lattice.setter
    def lattice(self, lattice):
        if lattice == 'strip' or lattice == '110' or lattice == 110:
            self._orient = 1
        elif lattice == 'quasi-mosaic' or lattice == '111' or lattice == 111:
            self._orient = 2
        else:
            raise ValueError(f"Illegal value {lattice} for 'lattice'! "
                            + "Only use 'strip' (110) or 'quasi-mosaic' (111).")
