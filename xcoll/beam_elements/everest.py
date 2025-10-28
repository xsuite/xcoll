# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt

from .base import BaseBlock, BaseCollimator, BaseCrystal, InvalidXcoll
from ..scattering_routines.everest import EverestEngine
from ..materials import Material, CrystalMaterial, db as material_db, _DEFAULT_MATERIAL, _DEFAULT_CRYSTALMATERIAL
from ..general import _pkg_root


# TODO:
#      We want these elements to behave as if 'iscollective = True' when doing twiss etc (because they would ruin the CO),
#      but as if 'iscollective = False' for normal tracking as it is natively in C...
#      Currently this is achieved with the hack '_tracking' which defaults to False after installation in the line, and is
#      only activated around the track command. Furthermore, because of 'iscollective = False' we need to specify
#      get_backtrack_element. We want it nicer..


class EverestBlock(BaseBlock):
    _xofields = {**BaseBlock._xofields,
        '_material':        Material,
        'rutherford_rng':   xt.RandomRutherford,
        '_tracking':        xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    allow_double_sided = False
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = {*BaseBlock._noexpr_fields, 'material'}
    _skip_in_to_dict       = [*BaseBlock._skip_in_to_dict, '_material']
    _store_in_to_dict      = [*BaseBlock._store_in_to_dict, 'material']
    _internal_record_class = BaseBlock._internal_record_class

    _depends_on = [BaseBlock, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','everest_block.h')
    ]

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
            kwargs['_material'] = _DEFAULT_MATERIAL
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


    @property
    def material(self):
        if self._material != _DEFAULT_MATERIAL:
            return self._material

    @material.setter
    def material(self, material):
        if material is None:
            material = _DEFAULT_MATERIAL
        elif isinstance(material, dict):
            material = Material.from_dict(material)
        elif isinstance(material, str):
            material = material_db[material]
        if not isinstance(material, Material) \
        or isinstance(material, CrystalMaterial):
            raise ValueError(f"Invalid material of type {type(material)}!")
        if self.material != material:
            self._material = material
            self.EverestBlock_set_material(el=self)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                                 _buffer=_buffer, _offset=_offset)


class EverestCollimator(BaseCollimator):
    _xofields = {**BaseCollimator._xofields,
        '_material':        Material,
        'rutherford_rng':   xt.RandomRutherford,
        '_tracking':        xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    allow_double_sided = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = {*BaseCollimator._noexpr_fields, 'material'}
    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_material']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material']
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, EverestEngine]

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
            kwargs['_material'] = _DEFAULT_MATERIAL
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)

    @property
    def material(self):
        if self._material != _DEFAULT_MATERIAL:
            return self._material

    @material.setter
    def material(self, material):
        if material is None:
            material = _DEFAULT_MATERIAL
        elif isinstance(material, dict):
            material = Material.from_dict(material)
        elif isinstance(material, str):
            material = material_db[material]
        if not isinstance(material, Material) \
        or isinstance(material, CrystalMaterial):
            raise ValueError(f"Invalid material of type {type(material)}!")
        if self.material != material:
            self._material = material
            self.EverestCollimator_set_material(el=self)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                            _buffer=_buffer, _offset=_offset)


class EverestCrystal(BaseCrystal):
    _xofields = {**BaseCrystal._xofields,
        'miscut':             xo.Float64,
        '_orient':            xo.Int8,
        '_critical_angle':    xo.Float64,
        '_critical_radius':   xo.Float64,
        '_material':          CrystalMaterial,
        'rutherford_rng':     xt.RandomRutherford,
        '_tracking':          xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    allow_double_sided = False
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = {*BaseCrystal._noexpr_fields, 'material', 'lattice'}
    _skip_in_to_dict       = [*BaseCrystal._skip_in_to_dict, '_orient', '_material']
    _store_in_to_dict      = [*BaseCrystal._store_in_to_dict, 'lattice', 'material']
    _internal_record_class = BaseCrystal._internal_record_class

    _depends_on = [BaseCrystal, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','everest_crystal.h')
    ]

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
            kwargs['_material'] = _DEFAULT_CRYSTALMATERIAL
            to_assign['lattice'] = kwargs.pop('lattice', 'strip')
            kwargs.setdefault('miscut', 0)
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


    @property
    def material(self):
        if self._material != _DEFAULT_CRYSTALMATERIAL:
            return self._material

    @material.setter
    def material(self, material):
        if material is None:
            material = _DEFAULT_CRYSTALMATERIAL
        elif isinstance(material, dict):
            material = Material.from_dict(material)
        elif isinstance(material, str):
            material = material_db[material]
        if not isinstance(material, CrystalMaterial):
            if isinstance(material, Material):
                if material.name == 'CarbonFibreCarbon':
                    material = material_db['CarbonCrystal']
                else:
                    material = material_db[f'{material.name}Crystal']
            else:
                raise ValueError(f"Invalid material of type {type(material)}!")
        if self.material != material:
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


    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                            _buffer=_buffer, _offset=_offset)


