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
#      We want these elements to behave as if 'iscollective = True' when doing twiss etc (because they would ruin the CO),
#      but as if 'iscollective = False' for normal tracking as it is natively in C...
#      Currently this is achieved with the hack '_tracking' which defaults to False after installation in the line, and is
#      only activated around the track command. Furthermore, because of 'iscollective = False' we need to specify 
#      get_backtrack_element. We want it nicer..

class EverestCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        '_material':        Material,
        'rutherford_rng':   xt.RandomRutherford,
        '_tracking':        xo.Int8
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [ *BaseCollimator._skip_in_to_dict, '_material' ]
    _store_in_to_dict      = [ *BaseCollimator._store_in_to_dict, 'material' ]
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','everest_collimator.h')
    ]


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
        super().__init__(**kwargs)
#         if '_xobject' not in kwargs:
#             self.random_generator.set_rutherford_by_xcoll_material(self.material)

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if not isinstance(material, Material):
            if not isinstance('material', dict) or material['__class__'] != "Material":
                raise ValueError("Invalid material!")
        self._material = material
#         self.random_generator.set_rutherford_by_xcoll_material(material)

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
        'miscut':         xo.Float64,
        '_orient':        xo.Int8,
        '_material':      CrystalMaterial,
        'rutherford_rng': xt.RandomRutherford,
        '_tracking':      xo.Int8
    }

    isthick = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_orient', '_material']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'lattice', 'material']
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','collimators_src','everest_crystal.h')
    ]


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            if kwargs.get('material') is None:
                raise ValueError("Need to provide a material to the collimator!")
            if not isinstance(kwargs['material'], CrystalMaterial):
                if not isinstance(kwargs['material'], dict) \
                or kwargs['material']['__class__'] != "CrystalMaterial":
                    raise ValueError("Invalid material!")
            kwargs['_material'] = kwargs.pop('material')
            kwargs.setdefault('bend', 0)
            kwargs.setdefault('xdim', 0)
            kwargs.setdefault('ydim', 0)
            kwargs.setdefault('thick', 0)
            kwargs.setdefault('miscut', 0)
            kwargs['_orient'] = _lattice_setter(kwargs.pop('lattice', 'strip'))
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
#         if '_xobject' not in kwargs:
#             self.random_generator.set_rutherford_by_xcoll_material(self.material)

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
        self._material = material
#         self.random_generator.set_rutherford_by_xcoll_material(material)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        # TODO: this should be an InvalidCollimator
        return xt.Drift(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)


def _lattice_setter(lattice):
    if lattice == 'strip' or lattice == '110' or lattice == 110:
        return 1
    elif lattice == 'quasi-mosaic' or lattice == '111' or lattice == 111:
        return 2
    else:
        raise ValueError(f"Illegal value {lattice} for 'lattice'! "
                        + "Only use 'strip' (110) or 'quasi-mosaic' (111).")
