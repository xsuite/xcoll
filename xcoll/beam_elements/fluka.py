# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from contextlib import contextmanager

import xobjects as xo
import xtrack as xt

from .base import BaseCollimator
from ..scattering_routines.fluka import track, FlukaEngine, assemblies, \
                    FlukaPrototype, FlukaAssembly, FlukaGenericAssembly, \
                    FlukaGenericCrystalAssembly
from ..scattering_routines.fluka.prototype import assemblies_wrong_jaw


class FlukaCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'fluka_id':              xo.Int16,    # Do not change! Should be 16 bit because of FlukaIO type
        'length_front':          xo.Float64,
        'length_back':           xo.Float64,
        '_tracking':             xo.Int8,
        '_acc_ionisation_loss':  xo.Float64,  # TODO: this is not very robust, for when a track is done with new particles etc
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = [ *BaseCollimator._store_in_to_dict, 'assembly' ]
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, FlukaEngine]

    _allowed_fields_when_frozen = ['_tracking', '_acc_ionisation_loss']

    def __new__(cls, *args, **kwargs):
        with cls._in_constructor():
            self = super().__new__(cls, *args, **kwargs)
        return self

    def __init__(self, **kwargs):
        with self.__class__._in_constructor():
            to_assign = {}
            if '_xobject' not in kwargs:
                kwargs.setdefault('_tracking', True)
                kwargs.setdefault('_acc_ionisation_loss', -1.)
                to_assign['name'] = FlukaEngine()._get_new_element_name()
                to_assign['assembly'] = kwargs.pop('assembly', None)
                if 'material' in kwargs or 'side' in kwargs:
                    if to_assign['assembly'] is not None:
                        raise ValueError('Cannot set both material/side and assembly!')
                    material = kwargs.pop('material', None)
                    if material is None:
                        raise ValueError('Need to provide material!')
                    length = kwargs.get('length', None)
                    if length is None:
                        raise ValueError('Need to provide length!')
                    side = kwargs.pop('side', 'both')
                    to_assign['assembly'] = FlukaGenericAssembly(material=material, side=side,
                                                                 length=length)
            super().__init__(**kwargs)
            for key, val in to_assign.items():
                setattr(self, key, val)
            if not hasattr(self, '_equivalent_drift'):
                self._equivalent_drift = xt.Drift(length=self.length)

    def __del__(self):
        self.assembly.remove_element(self.name, force=False)
        try:
            super().__del__()
        except AttributeError:
            pass

    def copy(self, **kwargs):
        obj = super().copy(**kwargs)
        obj.assembly = self.assembly
        return obj

    @property
    def material(self):
        if self.assembly is not None:
            return self.assembly.material

    @material.setter
    def material(self, material):
        self.assembly = FlukaGenericAssembly(material=material, side=self.side,
                                             length=self.length)

    @property
    def side(self):
        if self.assembly is not None:
            return self.assembly.side

    @side.setter
    def side(self, side):
        self.assembly = FlukaGenericAssembly(material=self.material, side=side,
                                             length=self.length)

    @property
    def assembly(self):
        for prototype in FlukaPrototype._assigned_registry.values():
            if self in prototype.elements:
                return prototype
        for prototype in FlukaAssembly._assigned_registry.values():
            if self in prototype.elements:
                return prototype
        return None

    @assembly.setter
    def assembly(self, assembly):
        if isinstance(assembly, str):
            if assembly in assemblies:
                assembly = assemblies[assembly]
            elif assembly in assemblies_wrong_jaw:
                print(f"Warning: Assembly '{assembly}' might be wrong.")
                assembly = assemblies_wrong_jaw[assembly]
            else:
                raise ValueError(f"Assembly (or prototype) '{assembly}' not present "
                               + f"in internal database. Please define it yourself.")
        elif not isinstance(assembly, FlukaPrototype) and assembly is not None:
            raise ValueError('Invalid assembly or prototype!')
        # Remove the element from the old assembly and add it to the new one
        if self.assembly:
            self.assembly.remove_element(self, force=False)
        if assembly:
            assembly.add_element(self, force=False)
            BaseCollimator.side.fset(self, self.assembly.side)
            if isinstance(self, BaseCrystal):
                BaseCrystal.bending_radius.fset(self, self.assembly.bending_radius)
        if self.assembly:
            if self.assembly.side is not None and self.assembly.side != self.side:
                self.side = self.assembly.side
                print(f"Warning: Side of collimator '{self.name}' was changed to '{self.side}' "
                    + f"to match the assembly '{self.assembly.name}'.")
            if self.assembly.length is not None:
                self.length_front = (self.assembly.length - self.length) / 2
                self.length_back = self.assembly.length - self.length - self.length_front

    def track(self, part):
        track(self, part)

    def __setattr__(self, name, value):
        if name not in self._allowed_fields_when_frozen \
        and FlukaEngine.is_running() is True:
            raise ValueError('Engine is running; FlukaCollimator is frozen.')
        super().__setattr__(name, value)

    @classmethod
    @contextmanager
    def _in_constructor(cls):
        original_setattr = cls.__setattr__
        def new_setattr(self, *args, **kwargs):
            return super().__setattr__( *args, **kwargs)
        cls.__setattr__ = new_setattr
        try:
            yield
        finally:
            cls.__setattr__ = original_setattr



class FlukaCrystal(BaseCrystal):
    _xofields = { **BaseCrystal._xofields,
        'fluka_id':              xo.Int16,    # Do not change! Should be 16 bit because of FlukaIO type
        'length_front':          xo.Float64,
        'length_back':           xo.Float64,
        '_tracking':             xo.Int8,
        '_acc_ionisation_loss':  xo.Float64,  # TODO: this is not very robust, for when a track is done with new particles etc
    }

    isthick = True
    allow_track = True
    iscollective = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    skip_in_loss_location_refinement = True

    _skip_in_to_dict       = BaseCrystal._skip_in_to_dict
    _store_in_to_dict      = [ *BaseCrystal._store_in_to_dict, 'assembly' ]
    _internal_record_class = BaseCrystal._internal_record_class

    _depends_on = [BaseCrystal, FlukaEngine]

    _allowed_fields_when_frozen = ['_tracking', '_acc_ionisation_loss']

    def __new__(cls, *args, **kwargs):
        with cls._in_constructor():
            self = super().__new__(cls, *args, **kwargs)
        return self

    def __init__(self, **kwargs):
        with self.__class__._in_constructor():
            to_assign = {}
            if '_xobject' not in kwargs:
                kwargs.setdefault('_tracking', True)
                kwargs.setdefault('_acc_ionisation_loss', -1.)
                to_assign['name'] = FlukaEngine()._get_new_element_name()
                to_assign['assembly'] = kwargs.pop('assembly', None)
                if to_assign['assembly']:
                    raise NotImplementedError('FlukaCrystalAssemblies not yet implemented!')
                material = None
                if 'material' in kwargs or 'side' in kwargs or 'bending_radius' in kwargs \
                or 'bending_angle' in kwargs:
                    if to_assign['assembly'] is not None:
                        raise ValueError('Cannot set both material and assembly!')
                    material = kwargs.pop('material', None)
                    if material is None:
                        raise ValueError('Need to provide material!')
                    length = kwargs.get('length', None)
                    if length is None:
                        raise ValueError('Need to provide length!')
                    bending_radius = kwargs.pop('bending_radius', None)
                    bending_angle = kwargs.pop('bending_angle', None)
                    if bending_radius is None and bending_angle is None:
                        raise ValueError('Need to provide bending radius or angle!')
                    if bending_radius is not None and bending_angle is not None:
                        raise ValueError('Cannot provide both bending radius and angle!')
                    side = kwargs.pop('side', None)
                    if side is None:
                        raise ValueError('Need to provide side!')
            super().__init__(**kwargs)
            for key, val in to_assign.items():
                setattr(self, key, val)
            if self.assembly is None and material:
                if bending_radius:
                    BaseCrystal.bending_radius.fset(self, bending_radius)
                if bending_angle:
                    BaseCrystal.bending_angle.fset(self, bending_angle)
                self.assembly = FlukaGenericCrystalAssembly(material=material, side=side,
                                         length=self.length, bending_radius=self.bending_radius)
            if not hasattr(self, '_equivalent_drift'):
                self._equivalent_drift = xt.Drift(length=self.length)

    def __del__(self):
        return FlukaCollimator.__del__(self, **kwargs)

    def copy(self, **kwargs):
        return FlukaCollimator.copy(self, **kwargs)

    @property
    def material(self):
        return FlukaCollimator.material.fget(self)

    @material.setter
    def material(self, material):
        self.assembly = FlukaGenericCrystalAssembly(material=material, side=self.side,
                               length=self.length, bending_radius=self.bending_radius)

    @property
    def side(self):
        return FlukaCollimator.side.fget(self)

    @side.setter
    def side(self, side):
        self.assembly = FlukaGenericCrystalAssembly(material=self.material, side=side,
                               length=self.length, bending_radius=self.bending_radius)

    @property
    def bending_radius(self):
        return BaseCrystal.bending_radius.fget(self)

    @bending_radius.setter
    def bending_radius(self, material):
        self.assembly = FlukaGenericCrystalAssembly(material=self.material, side=self.side,
                                         length=self.length, bending_radius=bending_radius)

    @property
    def assembly(self):
        return FlukaCollimator.assembly.fget(self)

    @assembly.setter
    def assembly(self, assembly):
        FlukaCollimator.assembly.fset(self, assembly)

    def track(self, part):
        FlukaCollimator.track(self, part)

    def __setattr__(self, name, value):
        if name not in self._allowed_fields_when_frozen \
        and FlukaEngine.is_running() is True:
            raise ValueError('Engine is running; FlukaCrystal is frozen.')
        super().__setattr__(name, value)

    @classmethod
    @contextmanager
    def _in_constructor(cls):
        original_setattr = cls.__setattr__
        def new_setattr(self, *args, **kwargs):
            return super().__setattr__( *args, **kwargs)
        cls.__setattr__ = new_setattr
        try:
            yield
        finally:
            cls.__setattr__ = original_setattr
