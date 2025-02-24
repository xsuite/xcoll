# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from contextlib import contextmanager

import xobjects as xo
import xtrack as xt

from .base import BaseCollimator
from ..scattering_routines.fluka import track, FlukaEngine, assemblies, FlukaPrototype
from ..scattering_routines.fluka.prototypes import assemblies_wrong_jaw


class FlukaCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'fluka_id':              xo.Int16,    # Do not change! Should be 16 bit because of FlukaIO type
        'length_front':          xo.Float64,
        'length_back':           xo.Float64,
        '_tracking':             xo.Int8,
        '_acc_ionisation_loss':  xo.Float64,  # TODO: this is not very robust, for when a track is done with new particles etc
        '_assembly':             FlukaPrototype
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
                to_assign['assembly'] = kwargs.pop('assembly', None)
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
    def name(self):
        return BaseCollimator.name.fget(self)

    @name.setter
    def name(self, val):
        if self.name is not None:
            self.assembly.remove_element(self.name, force=False)
        BaseCollimator.name.fset(self, val)
        if self.name is not None:
            self.assembly.add_element(self.name, force=False)

    @property
    def assembly(self):
        if not hasattr(self, '_assembly'):
            self._assembly = FlukaPrototype()
        return self._assembly

    @assembly.setter
    def assembly(self, assembly):
        if assembly is None:
            assembly = FlukaPrototype()
        if isinstance(assembly, str):
            if assembly in assemblies:
                assembly = assemblies[assembly]
            elif assembly in assemblies_wrong_jaw:
                print(f"Warning: Assembly '{assembly}' might be wrong.")
                assembly = assemblies_wrong_jaw[assembly]
            else:
                raise ValueError(f"Assembly (or prototype) '{assembly}' not present "
                               + f"in internal database. Please define it yourself.")
        elif not isinstance(assembly, FlukaPrototype):
            raise ValueError('Invalid assembly or prototype!')
        if self.name is not None:
            # Remove the element from the old assembly and add it to the new one
            self.assembly.remove_element(self.name, force=False)
            assembly.add_element(self.name, force=False)
        self._assembly = assembly
        if self.assembly.side is not None and self.assembly.side != self.side:
            self.side = self.assembly.side
            print(f"Warning: Side of collimator '{self.name}' was changed to '{self.side}' "
                + f"to match the assembly '{self.assembly.name}'.")

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
