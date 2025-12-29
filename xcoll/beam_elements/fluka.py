# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from contextlib import contextmanager

import xobjects as xo
import xtrack as xt

from .base import BaseCollimator, BaseCrystal
from ..general import _pkg_root
from ..scattering_routines.fluka import track_pre, track_core, track_post, FlukaEngine, \
                                        FlukaPrototype, FlukaAssembly, create_generic_assembly
from ..materials import _resolve_material
from ..constants import HIT_ON_FLUKA_COLL


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

    _depends_on = [BaseCollimator, FlukaEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','fluka_collimator.h')
    ]

    _noexpr_fields         = {*BaseCollimator._noexpr_fields, 'material', 'assembly'}
    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_tracking', '_acc_ionisation_loss']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material', 'assembly', 'height', 'width', 'side']
    _internal_record_class = BaseCollimator._internal_record_class

    _allowed_fields_when_frozen = ['_tracking', '_acc_ionisation_loss']

    def __new__(cls, *args, **kwargs):
        with cls._in_constructor():
            self = super().__new__(cls, *args, **kwargs)
        return self

    def __init__(self, **kwargs):
        import xcoll as xc
        with self.__class__._in_constructor(self):
            to_assign = {}
            if '_xobject' not in kwargs:
                kwargs.setdefault('_tracking', True)
                kwargs.setdefault('_acc_ionisation_loss', -1.)
                to_assign['name'] = xc.fluka.engine._get_new_element_name()
                assembly = kwargs.pop('assembly', None)
                material = _resolve_material(kwargs.pop('material', None), ref='fluka', allow_none=True)
                side = kwargs.pop('side', None)
                width = kwargs.pop('width', None)
                height = kwargs.pop('height', None)
                if assembly is not None:
                    # Use the provided assembly, check consistency later
                    generic = False
                    to_assign['assembly'] = assembly
                else:
                    # Create a generic assembly
                    generic = True
                    if material is None:
                        raise ValueError('Need to provide material!')
                    length = kwargs.get('length', None)
                    if length is None:
                        raise ValueError('Need to provide length!')

            super().__init__(**kwargs)
            for key, val in to_assign.items():
                setattr(self, key, val)
            if generic:
                side = self._get_side_from_input(side)
                self.assembly = create_generic_assembly(material=material,
                                    side=side, length=self.length, width=width,
                                    height=height)
            else:
                # Check consistency
                if self.assembly.material is not None and material is not None \
                and self.assembly.material != material:
                    print("Material in assembly and provided `material` "
                        + "argument do not agree! Ignored the latter.")
                if self.assembly.side is not None and side is not None \
                and self.assembly.side != side:
                    print("Side in assembly and provided `side` "
                        + "argument do not agree! Ignored the latter.")
                if self.assembly.width is not None and width is not None \
                and self.assembly.width != width:
                    print("Width in assembly and provided `width` "
                        + "argument do not agree! Ignored the latter.")
                if self.assembly.height is not None and height is not None \
                and self.assembly.height != height:
                    print("Height in assembly and provided `height` "
                        + "argument do not agree! Ignored the latter.")
            if not hasattr(self, '_equivalent_drift'):
                self._equivalent_drift = xt.Drift(length=self.length)

    def __del__(self):
        if self.assembly:
            self.assembly.remove_element(self, force=False)
        try:
            super().__del__()
        except (AttributeError, TypeError): # During shutdown, if the parent class is already garbage collected
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
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change material of non-generic assembly!')
            material = _resolve_material(material, ref='fluka', allow_none=False)
            if self.material != material:
                self.assembly = create_generic_assembly(material=material,
                                side=self.side, length=self.length, width=self.width,
                                height=self.height)

    @property
    def height(self):
        if self.assembly is not None:
            return self.assembly.height

    @height.setter
    def height(self, height):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change height of non-generic assembly!')
            self.assembly = create_generic_assembly(material=self.material,
                            side=self.side, length=self.length, width=self.width,
                            height=height)

    @property
    def width(self):
        if self.assembly is not None:
            return self.assembly.width

    @width.setter
    def width(self, width):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change width of non-generic assembly!')
            self.assembly = create_generic_assembly(material=self.material,
                            side=self.side, length=self.length, width=width,
                            height=self.height)
    @property
    def side(self):
        if self.assembly is not None:
            return self.assembly.side

    @side.setter
    def side(self, side):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change side of non-generic assembly!')
            side = self._get_side_from_input(side)
            self.assembly = create_generic_assembly(material=self.material,
                            side=side, length=self.length, width=self.width,
                            height=self.height)

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
    def assembly(self, val):
        import xcoll as xc
        if isinstance(val, str):
            if val in xc.fluka.assemblies:
                val = xc.fluka.assemblies[val]
            elif val in xc.fluka.prototypes:
                val = xc.fluka.prototypes[val]
            else:
                raise ValueError(f"Unknown assembly or prototype '{val}'.")
            if val.is_broken:
                print(f'Warning: assembly or prototype {val.name} is broken!')
        elif not isinstance(val, FlukaPrototype) and val is not None:
            raise ValueError(f'Invalid assembly or prototype {val}!')
        # Remove the element from the old assembly and add it to the new one
        if self.assembly:
            self.assembly.remove_element(self, force=False)
        if val:
            val.add_element(self, force=False)
        if self.assembly:
            if self.assembly.length is not None:
                self.length_front = (self.assembly.length - self.length) / 2
                self.length_back = self.assembly.length - self.length - self.length_front
            if self.assembly.side is not None:
                self._get_side_from_input(self.assembly.side)


    def enable_scattering(self):
        import xcoll as xc
        xc.fluka.environment.assert_environment_ready()
        if not xc.fluka.engine.is_running():
            raise RuntimeError("FLUKA engine is not running.")
        super().enable_scattering()

    def track(self, part):
        if track_pre(self, part):
            if self.material != "vacuum":
                super().track(part)
            else:
                part.state[part.state == 1] = HIT_ON_FLUKA_COLL
            track_core(self, part)
            if self.material != "vacuum":
                track_post(self, part)
        else:
            self._drift(part)

    def _drift(self, particles, length=None):
        if length is None:
            length = self.length
        if length != self.length:
            old_length = self._equivalent_drift.length
            self._equivalent_drift.length = length
        self._equivalent_drift.track(particles)
        if length != self.length:
            self._equivalent_drift.length = old_length

    def __setattr__(self, name, value):
        import xcoll as xc
        if name not in self._allowed_fields_when_frozen \
        and xc.fluka.engine.is_running() is True:
            raise ValueError('Engine is running; FlukaCollimator is frozen.')
        super().__setattr__(name, value)

    @classmethod
    @contextmanager
    def _in_constructor(cls, self=None):
        original_setattr = cls.__setattr__
        if self is not None:
            self._being_constructed_ = True
        def new_setattr(self, *args, **kwargs):
            return super().__setattr__( *args, **kwargs)
        cls.__setattr__ = new_setattr
        try:
            yield
        finally:
            cls.__setattr__ = original_setattr
            if self is not None:
                self._being_constructed_ = False

    def _being_constructed(self):
        if hasattr(self, '_being_constructed_'):
            return self._being_constructed_
        else:
            return False

    # ===================================================
    # ===   Hacks to use parent setters and getters   ===
    # ===================================================

    def _get_side_from_input(self, side):
        # Set / get using _side field (even though it is not used directly)
        if side is not None:
            BaseCollimator.side.fset(self, side)
            side = BaseCollimator.side.fget(self)
            return side



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

    _noexpr_fields         = {*BaseCrystal._noexpr_fields, 'material', 'assembly'}
    _skip_in_to_dict       = [*BaseCrystal._skip_in_to_dict, '_tracking', '_acc_ionisation_loss']
    _store_in_to_dict      = [*BaseCrystal._store_in_to_dict, 'material', 'assembly', 'height', 'width', 'side']
    _internal_record_class = BaseCrystal._internal_record_class

    _depends_on = [BaseCrystal, FlukaEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','fluka_crystal.h')
    ]

    _allowed_fields_when_frozen = ['_tracking', '_acc_ionisation_loss']

    def __new__(cls, *args, **kwargs):
        with cls._in_constructor():
            self = super().__new__(cls, *args, **kwargs)
        return self

    def __init__(self, **kwargs):
        with self.__class__._in_constructor(self):
            import xcoll as xc
            to_assign = {}
            generic = False
            if '_xobject' not in kwargs:
                kwargs.setdefault('_tracking', True)
                kwargs.setdefault('_acc_ionisation_loss', -1.)
                to_assign['name'] = xc.fluka.engine._get_new_element_name()
                assembly = kwargs.pop('assembly', None)
                if assembly:
                    raise NotImplementedError('FlukaCrystalAssemblies not yet implemented!')
                # kwargs for generic assembly creation
                material = _resolve_material(kwargs.pop('material', None), ref='fluka', allow_none=True)
                side = kwargs.pop('side', None)
                width = kwargs.pop('width', None)
                height = kwargs.pop('height', None)
                bending_radius = kwargs.pop('bending_radius', None)
                bending_angle = kwargs.pop('bending_angle', None)
                if assembly is not None:
                    # Use the provided assembly, check consistency later
                    generic = False
                    to_assign['assembly'] = assembly
                else:
                    # Create a generic assembly
                    generic = True
                    if material is None:
                        raise ValueError('Need to provide material!')
                    length = kwargs.get('length', None)
                    if length is None:
                        raise ValueError('Need to provide length!')
                    if side is None:
                        raise ValueError('Need to provide side!')
                    if bending_radius is None and bending_angle is None:
                        raise ValueError('Need to provide bending radius or angle!')
                    elif bending_radius is not None and bending_angle is not None:
                        raise ValueError('Cannot provide both bending radius and angle!')
                # Some defaults to keep BaseCrystal happy
                kwargs.setdefault('_bending_radius', 1.)
                length = kwargs.get('length', kwargs['length'] or 1 )
                kwargs.setdefault('_bending_angle', np.arcsin(length))
            super().__init__(**kwargs)
            for key, val in to_assign.items():
                setattr(self, key, val)
            if generic:
                if bending_radius is None:
                    bending_radius = self._get_bending_radius_from_angle(bending_angle)
                else:
                    self._get_bending_angle_from_radius(bending_radius)  # To set internal fields correctly
                side = self._get_side_from_input(side)
                self.assembly = create_generic_assembly(is_crystal=True, material=material,
                                side=side, length=self.length, width=width, height=height,
                                bending_radius=bending_radius)
            else:
                # Check consistency
                if self.assembly.material is not None and material is not None \
                and self.assembly.material != material:
                    print("Material in assembly and provided `material` "
                        + "argument do not agree! Ignored the latter.")
                if self.assembly.side is not None and side is not None \
                and self.assembly.side != side:
                    print("Side in assembly and provided `side` "
                        + "argument do not agree! Ignored the latter.")
                if self.assembly.width is not None and width is not None \
                and self.assembly.width != width:
                    print("Width in assembly and provided `width` "
                        + "argument do not agree! Ignored the latter.")
                if self.assembly.height is not None and height is not None \
                and self.assembly.height != height:
                    print("Height in assembly and provided `height` "
                        + "argument do not agree! Ignored the latter.")
                if self.assembly.bending_radius is not None and bending_radius is not None \
                and self.assembly.bending_radius != bending_radius:
                    print("Bending_radius in assembly and provided `bending_radius` "
                        + "argument do not agree! Ignored the latter.")
                # TODO: check bending_angle consistency too
            if not hasattr(self, '_equivalent_drift'):
                self._equivalent_drift = xt.Drift(length=self.length)

    def __del__(self, **kwargs):
        FlukaCollimator.__del__(self, **kwargs)

    def copy(self, **kwargs):
        return FlukaCollimator.copy(self, **kwargs)

    @property
    def material(self):
        return FlukaCollimator.material.fget(self)

    @material.setter
    def material(self, material):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change material of non-generic assembly!')
            material = _resolve_material(material, ref='fluka', allow_none=False)
            if self.material != material:
                self.assembly = create_generic_assembly(is_crystal=True, material=material,
                                side=self.side, length=self.length, width=self.width,
                                height=self.height, bending_radius=self.bending_radius)

    @property
    def height(self):
        return FlukaCollimator.height.fget(self)

    @height.setter
    def height(self, height):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change height of non-generic assembly!')
            self.assembly = create_generic_assembly(is_crystal=True, material=self.material,
                            side=self.side, length=self.length, width=self.width,
                            height=height, bending_radius=self.bending_radius)

    @property
    def width(self):
        return FlukaCollimator.width.fget(self)

    @width.setter
    def width(self, width):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change width of non-generic assembly!')
            self.assembly = create_generic_assembly(is_crystal=True, material=self.material,
                            side=self.side, length=self.length, width=width,
                            height=self.height, bending_radius=self.bending_radius)

    @property
    def side(self):
        if self.assembly is not None:
            return FlukaCollimator.side.fget(self)

    @side.setter
    def side(self, side):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change side of non-generic assembly!')
            side = self._get_side_from_input(side)
            self.assembly = create_generic_assembly(is_crystal=True, material=self.material,
                            side=side, length=self.length, width=self.width,
                            height=self.height, bending_radius=self.bending_radius)
    @property
    def bending_radius(self):
        if self.assembly is not None:
            return self.assembly.bending_radius

    @bending_radius.setter
    def bending_radius(self, bending_radius):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change bending radius of non-generic assembly!')
            self._get_bending_angle_from_radius(bending_radius) # To set internal fields correctly
            self.assembly = create_generic_assembly(is_crystal=True, material=self.material,
                            side=self.side, length=self.length, width=self.width,
                            height=self.height, bending_radius=bending_radius)

    @property
    def bending_angle(self):
        return self._get_bending_angle_from_radius(self.bending_radius)

    @bending_angle.setter
    def bending_angle(self, bending_angle):
        if not self._being_constructed():
            if self.assembly.fedb_series != 'generic':
                raise ValueError('Cannot change bending angle of non-generic assembly!')
            self.bending_radius = self._get_bending_radius_from_angle(bending_angle)

    @property
    def assembly(self):
        return FlukaCollimator.assembly.fget(self)

    @assembly.setter
    def assembly(self, assembly):
        FlukaCollimator.assembly.fset(self, assembly)


    def enable_scattering(self):
        import xcoll as xc
        xc.fluka.environment.assert_environment_ready()
        if not xc.fluka.engine.is_running():
            raise RuntimeError("FLUKA engine is not running.")
        super().enable_scattering()

    def track(self, part):
        if track_pre(self, part):
            # super().track(part)
            part.state[part.state == 1] = HIT_ON_FLUKA_COLL
            track_core(self, part)
            track_post(self, part)
        else:
            self._drift(part)

    def _drift(self, particles, length=None):
        if length is None:
            length = self.length
        if length != self.length:
            old_length = self._equivalent_drift.length
            self._equivalent_drift.length = length
        self._equivalent_drift.track(particles)
        if length != self.length:
            self._equivalent_drift.length = old_length

    def __setattr__(self, name, value):
        import xcoll as xc
        if name not in self._allowed_fields_when_frozen \
        and xc.fluka.engine.is_running() is True:
            raise ValueError('Engine is running; FlukaCrystal is frozen.')
        super().__setattr__(name, value)

    _in_constructor = FlukaCollimator._in_constructor
    _being_constructed = FlukaCollimator._being_constructed


    # ===================================================
    # ===   Hacks to use parent setters and getters   ===
    # ===================================================

    def _get_side_from_input(self, side):
        # Set / get using _side field (even though it is not used directly)
        if side is not None:
            BaseCrystal.side.fset(self, side)
            side = BaseCrystal.side.fget(self)
            return side

    def _get_bending_radius_from_angle(self, bending_angle):
        # Set / get using _bending_radius and _bending_angle fields (even though they are not used directly)
        if self.assembly:
            old_length = self.length
            self.length = self.assembly.length
        BaseCrystal.bending_angle.fset(self, bending_angle)
        bending_radius = BaseCrystal.bending_radius.fget(self)
        if self.assembly:
            self.length = old_length
        return bending_radius

    def _get_bending_angle_from_radius(self, bending_radius):
        # Set / get using _bending_radius and _bending_angle fields (even though they are not used directly)
        if self.assembly:
            old_length = self.length
            self.length = self.assembly.length
        BaseCrystal.bending_radius.fset(self, bending_radius)
        bending_angle = BaseCrystal.bending_angle.fget(self)
        if self.assembly:
            self.length = old_length
        return bending_angle
