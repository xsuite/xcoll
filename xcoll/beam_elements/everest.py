# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xpart as xp
import xtrack as xt

from .base import BaseBlock, BaseCollimator, InvalidXcoll
from ..scattering_routines.everest import GeneralMaterial, Material, CrystalMaterial, EverestEngine
from ..general import _pkg_root


# TODO:
#      We want these elements to behave as if 'iscollective = True' when doing twiss etc (because they would ruin the CO),
#      but as if 'iscollective = False' for normal tracking as it is natively in C...
#      Currently this is achieved with the hack '_tracking' which defaults to False after installation in the line, and is
#      only activated around the track command. Furthermore, because of 'iscollective = False' we need to specify
#      get_backtrack_element. We want it nicer..

# TODO: _per_particle_kernels should be a normal kernel (such that we don't need to pass a dummy Particles() )

class EverestBlock(BaseBlock):
    _xofields = { **BaseBlock._xofields,
        '_material':        Material,
        'rutherford_rng':   xt.RandomRutherford,
        '_tracking':        xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = ['_material']
    _store_in_to_dict      = ['material']
    _internal_record_class = BaseBlock._internal_record_class

    _depends_on = [BaseBlock, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','everest_block.h')
    ]

    _kernels = {
        'EverestBlock_set_material': xo.Kernel(
                c_name='EverestBlock_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            mat = kwargs.pop('material', None)
            if mat is None:
                raise ValueError("Need to provide a material to the block!")
            if not isinstance(mat, Material):
                if not isinstance(mat, dict) \
                or mat['__class__'] != "Material":
                    raise ValueError("Invalid material!")
            kwargs['_material'] = mat
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
            use_prebuilt_kernels = kwargs.pop('use_prebuilt_kernels', True)
        super().__init__(**kwargs)
        if '_xobject' not in kwargs:
            # TODO: temporarily until release v0.3.4
            self.compile_kernels(only_if_needed=True)
            self._context.kernels.EverestBlock_set_material(el=self)


    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if not isinstance(material, Material):
            if not isinstance('material', dict) or material['__class__'] != "Material":
                raise ValueError("Invalid material!")
        if not xt.line._dicts_equal(self.material.to_dict(), material.to_dict()):
            self._material = material
            # TODO: temporarily until release v0.3.4
            self.compile_kernels(only_if_needed=True)
            self._context.kernels.EverestBlock_set_material(el=self)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                                 _buffer=_buffer, _offset=_offset)


class EverestCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        '_material':        Material,
        'rutherford_rng':   xt.RandomRutherford,
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

    _depends_on = [BaseCollimator, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','everest_collimator.h')
    ]

    _kernels = {
        'EverestCollimator_set_material': xo.Kernel(
                c_name='EverestCollimator_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            if kwargs.get('material') is None:
                raise ValueError("Need to provide a material to the collimator!")
            if not isinstance(kwargs['material'], Material):
                if not isinstance(kwargs['material'], dict) \
                or kwargs['material']['__class__'] != "Material":
                    raise ValueError("Invalid material!")
            kwargs['_material'] = kwargs.pop('material')
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
            use_prebuilt_kernels = kwargs.pop('use_prebuilt_kernels', True)
        super().__init__(**kwargs)
        if '_xobject' not in kwargs:
            # TODO: temporarily until release v0.3.4
            self.compile_kernels(only_if_needed=True)
            self._context.kernels.EverestCollimator_set_material(el=self)

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if not isinstance(material, Material):
            if not isinstance('material', dict) or material['__class__'] != "Material":
                raise ValueError("Invalid material!")
        if not xt.line._dicts_equal(self.material.to_dict(), material.to_dict()):
            self._material = material
            # TODO: temporarily until release v0.3.4
            self.compile_kernels(only_if_needed=True)
            self._context.kernels.EverestCollimator_set_material(el=self)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                                 _buffer=_buffer, _offset=_offset)



class EverestCrystal(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'align_angle':        xo.Float64,  #  = - sqrt(eps/beta)*alpha*nsigma
        '_bending_radius':    xo.Float64,
        '_bending_angle':     xo.Float64,
        '_critical_angle':    xo.Float64,
        'xdim':               xo.Float64,
        'ydim':               xo.Float64,
        'thick':              xo.Float64,
        'miscut':             xo.Float64,
        '_orient':            xo.Int8,
        '_material':          CrystalMaterial,
        'rutherford_rng':     xt.RandomRutherford,
        '_tracking':          xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_orient', '_material', '_bending_radius',
                              '_bending_angle']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'lattice', 'material', 'bending_radius', 'bending_angle']
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','everest_crystal.h')
    ]

    _kernels = {
        'EverestCrystal_set_material': xo.Kernel(
                c_name='EverestCrystal_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            if kwargs.get('material') is None:
                raise ValueError("Need to provide a material to the collimator!")
            if not isinstance(kwargs['material'], CrystalMaterial):
                if not isinstance(kwargs['material'], dict) \
                or kwargs['material']['__class__'] != "CrystalMaterial":
                    raise ValueError("Invalid material!")
            kwargs['_material'] = kwargs.pop('material')
            bending_radius = False
            bending_angle  = False
            if 'bending_radius' in kwargs:
                if 'bending_angle' in kwargs:
                    raise ValueError("Need to choose between 'bending_radius' and 'bending_angle'!")
                bending_radius = kwargs['bending_radius']
            elif 'bending_angle' in kwargs:
                bending_angle = kwargs['bending_angle']
            kwargs['_bending_radius'] = kwargs.pop('bending_radius',0)
            kwargs['_bending_angle'] = kwargs.pop('bending_angle', 0)
            kwargs.setdefault('xdim', 0)
            kwargs.setdefault('ydim', 0)
            kwargs.setdefault('thick', 0)
            kwargs.setdefault('miscut', 0)
            kwargs['_orient'] = _lattice_setter(kwargs.pop('lattice', 'strip'))
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
            use_prebuilt_kernels = kwargs.pop('use_prebuilt_kernels', True)
        super().__init__(**kwargs)
        if '_xobject' not in kwargs:
            if bending_radius:
                self._bending_angle = np.arcsin(self.active_length/bending_radius)
            if bending_angle:
                self._bending_radius = self.active_length / np.sin(bending_angle)
            # TODO: temporarily until release v0.3.4
            self.compile_kernels(only_if_needed=True)
            self._context.kernels.EverestCrystal_set_material(el=self)


    @property
    def critical_angle(self):
        return self._critical_angle if abs(self._critical_angle) > 1.e-10 else None

    @property
    def bending_radius(self):
        return self._bending_radius

    @bending_radius.setter
    def bending_radius(self, bending_radius):
        self._bending_radius = bending_radius
        self._bending_angle = np.arcsin(self.active_length/bending_radius)

    @property
    def bending_angle(self):
        return self._bending_angle

    @bending_angle.setter
    def bending_angle(self, bending_angle):
        self._bending_angle = bending_angle
        self._bending_radius = self.active_length / np.sin(bending_angle)

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
        self._orient = _lattice_setter(lattice)

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if not isinstance(material, CrystalMaterial):
            if not isinstance(material, dict) or material['__class__'] != "CrystalMaterial":
                raise ValueError("Invalid material!")
        if not xt.line._dicts_equal(self.material.to_dict(), material.to_dict()):
            self._material = material
            # TODO: temporarily until release v0.3.4
            self.compile_kernels(only_if_needed=True)
            self._context.kernels.EverestCrystal_set_material(el=self)


    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                                 _buffer=_buffer, _offset=_offset)


def _lattice_setter(lattice):
    if lattice == 'strip' or lattice == '110' or lattice == 110:
        return 1
    elif lattice == 'quasi-mosaic' or lattice == '111' or lattice == 111:
        return 2
    else:
        raise ValueError(f"Illegal value {lattice} for 'lattice'! "
                        + "Only use 'strip' (110) or 'quasi-mosaic' (111).")

