import numpy as np

import xobjects as xo

from xcoll import BaseCollimator, Material, CrystalMaterial
from ..scattering_routines.pyeverest import track


# TODO: remove dx, dy, offset, tilt, as this should only be in colldb (and here only the jaw positions)
class PyEverestCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        '_material':   Material,
        '_tracking':  xo.Int8
    }

    isthick = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class


    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            if kwargs.get('material') is None:
                raise ValueError("Need to provide a material to the collimator!")
            if not isinstance(kwargs['material'], Material):
                if not isinstance(kwargs['material'], dict) \
                or kwargs['material']['__class__'] != "Material":
                    raise ValueError("Invalid material!")
            kwargs['_material'] = kwargs.pop('material')
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)


    def track(self, particles):  # TODO: write impacts
        track(self, particles)
        return

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if not isinstance(material, Material):
            if not isinstance('material', dict) or material['__class__'] != "Material":
                raise ValueError("Invalid material!")
        self._material = material



class PyEverestCrystal(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'align_angle': xo.Float64,  #  = - sqrt(eps/beta)*alpha*nsigma
        'bend':        xo.Float64,
        'xdim':        xo.Float64,
        'ydim':        xo.Float64,
        'thick':       xo.Float64,
        'crytilt':     xo.Float64,
        'miscut':      xo.Float64,
        '_orient':      xo.Float64,
        '_material':    CrystalMaterial,
        '_tracking':   xo.Int8
    }

    isthick = True
    iscollective = True
    behaves_like_drift = True
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class


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
            kwargs.setdefault('crytilt', 0)
            kwargs.setdefault('miscut', 0)
            kwargs['_orient'] = _lattice_setter(kwargs.pop('lattice', 'strip'))
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)


    def track(self, particles):  # TODO: write impacts
        track(self, particles)
        return

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


def _lattice_setter(lattice):
    if lattice == 'strip' or lattice == '110' or lattice == 110:
        return 1
    elif lattice == 'quasi-mosaic' or lattice == '111' or lattice == 111:
        return 2
    else:
        raise ValueError(f"Illegal value {lattice} for 'lattice'! "
                        + "Only use 'strip' (110) or 'quasi-mosaic' (111).")

