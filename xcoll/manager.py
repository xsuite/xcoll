# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import json
from pathlib import Path
import numpy as np
import pandas as pd

from .beam_elements import BaseCollimator, BlackAbsorber, EverestCollimator, EverestCrystal, collimator_classes, element_classes
from .install import install_elements
from .colldb import CollimatorDatabase
from .impacts import CollimatorImpacts
from .scattering_routines.everest.materials import SixTrack_to_xcoll, CrystalMaterial

import xobjects as xo
import xpart as xp
import xtrack as xt


# Logic of buffers:
#    A buffer is a continuous stream of memory, containing a bunch of xobjects.
#    It is an architecture-invariant way to work with pointers, as the pointers to the individual
#    xobjects are stored as offsets to the start of the buffer.
#    The tracker manages three buffers: one for the line and all elements, one for the particles,
#    and an io_buffer for record logging.

class CollimatorManager:

    _init_vars = ['_colldb', 'line', 'beam', 'capacity', 'record_impacts', '_context', '_buffer', 'io_buffer', '_part_buffer']
    _init_var_defaults = {'_colldb': None, 'beam': None, 'capacity': 1e6, 'record_impacts': False, \
                          '_context': None, '_buffer': None, 'io_buffer': None, '_part_buffer': None}

    # -------------------------------
    # ------ Loading functions ------
    # -------------------------------

    @classmethod
    def from_yaml(cls, file, **kwargs):
        if '_colldb' in kwargs.keys():
            raise ValueError("Cannot load CollimatorDatabase file and specify '_colldb' argument in loader!")
        kwargs_colldb  = {key: val for key, val in kwargs.items() if key in CollimatorDatabase._init_vars}
        kwargs_manager = {key: val for key, val in kwargs.items() if key in cls._init_vars}
        colldb = CollimatorDatabase.from_yaml(file, **kwargs_colldb)
        return cls(_colldb=colldb, **kwargs_manager)

    @classmethod
    def from_json(cls, file, **kwargs):
        if '_colldb' in kwargs.keys():
            raise ValueError("Cannot load CollimatorDatabase file and specify '_colldb' argument in loader!")
        kwargs_colldb  = {key: val for key, val in kwargs.items() if key in CollimatorDatabase._init_vars}
        kwargs_manager = {key: val for key, val in kwargs.items() if key in cls._init_vars}
        colldb = CollimatorDatabase.from_json(file, **kwargs_colldb)
        return cls(_colldb=colldb, **kwargs_manager)

    @classmethod
    def from_dict(cls, file, **kwargs):
        if '_colldb' in kwargs.keys():
            raise ValueError("Cannot load CollimatorDatabase file and specify '_colldb' argument in loader!")
        kwargs_colldb  = {key: val for key, val in kwargs.items() if key in CollimatorDatabase._init_vars}
        kwargs_manager = {key: val for key, val in kwargs.items() if key in cls._init_vars}
        colldb = CollimatorDatabase.from_dict(file, **kwargs_colldb)
        return cls(_colldb=colldb, **kwargs_manager)

    @classmethod
    def from_SixTrack(cls, file, **kwargs):
        if '_colldb' in kwargs.keys():
            raise ValueError("Cannot load CollimatorDatabase file and specify '_colldb' argument in loader!")
        kwargs_colldb  = {key: val for key, val in kwargs.items() if key in CollimatorDatabase._init_vars}
        kwargs_manager = {key: val for key, val in kwargs.items() if key in cls._init_vars}
        colldb = CollimatorDatabase.from_SixTrack(file, **kwargs_colldb)
        return cls(_colldb=colldb, **kwargs_manager)


    def __init__(self, **kwargs):
        # Get all arguments
        for var in self._init_vars:
            if var in self._init_var_defaults:
                kwargs.setdefault(var, self._init_var_defaults[var])
            elif var not in kwargs.keys():
                raise ValueError(f"CollimatorManager is missing required argument '{var}'!")

        colldb = kwargs['_colldb']
        if not isinstance(colldb, CollimatorDatabase):
            # TODO allow None as empty database
            raise ValueError("The variable '_colldb' needs to be an xcoll CollimatorDatabase object!")
        else:
            self.colldb = colldb

        line = kwargs['line']
        beam = kwargs['beam']
        if isinstance(beam, str):
            beam = int(beam[-1])
        if not isinstance(line, xt.Line):
            raise ValueError("The variable 'line' needs to be an xtrack Line object!")
        else:
            self.line = line
        self.line._needs_rng = True  # TODO not needed if only BlackAbsorbers
        if beam is not None and beam > 1:
            self._line_is_reversed = True
        else:
            self._line_is_reversed = False
        self._machine_length = self.line.get_length()

        # Create _buffer, _context, and _io_buffer
        _buffer  = kwargs['_buffer']
        _context = kwargs['_context']
        if _buffer is None:
            if _context is None:
                _context = xo.ContextCpu()
            _buffer = _context.new_buffer()
        elif _context is not None and _buffer.context != _context:
            raise ValueError("The provided buffer and context do not match! "
                             + "Make sure the buffer is generated inside the provided context, or alternatively, "
                             + "only pass one of _buffer or _context.")
        self._buffer = _buffer
        self._part_buffer  = kwargs['_part_buffer']
        if self._part_buffer is not None and self._part_buffer.context != _context:
            raise ValueError("The provided particle buffer and context do not match! "
                             + "Make sure the particle buffer is generated inside the provided context.")

        # TODO: currently capacity is only for io_buffer (hence for impacts).
        # Do we need it in the _buffer as well?
        self._capacity = int(kwargs['capacity'])
        io_buffer = kwargs['io_buffer']
        if io_buffer is None:
            io_buffer = xt.new_io_buffer(_context=self._buffer.context, capacity=self.capacity)
        elif self._buffer.context != io_buffer._context:
            raise ValueError("The provided io_buffer lives on a different context than the buffer!")
        self._io_buffer = io_buffer

        # Initialise impacts table
        self._recording_elements = []
        self.record_impacts = kwargs['record_impacts']
        self._impacts = None

        # Variables for lossmap
        self._lossmap = None
        self._summary = None
        self._part    = None


    def __getitem__(self, name):
        return self.colldb[name]


    @property
    def machine_length(self):
        return self._machine_length

    @property
    def impacts(self):
        return self._impacts

    @property
    def record_impacts(self):
        return self._record_impacts

    @record_impacts.setter
    def record_impacts(self, record_impacts):
        if record_impacts is False or record_impacts is None:
            record_stop  = [self.line[name] for name in self._recording_elements]
            record_start = False
        elif record_impacts is True:
            record_stop  = []
            record_start = True
        else:
            if not hasattr(record_impacts, '__iter__') or isinstance(record_impacts, str):
                record_impacts = [record_impacts]
            record_stop  = [self.line[name] for name in set(self._recording_elements) - set(record_impacts)]
            record_start = record_impacts

        # We stop logging these elements
        if record_stop:
            if self.line.tracker is not None:
                self.line.tracker._check_invalidated()
            xt.stop_internal_logging(elements=record_stop)
            # Removed the stopped collimators from list of logged elements
            self._recording_elements = list(set(self._recording_elements) - set(record_stop))

        # We will start logging these elements
        self._record_impacts = record_start
        self._set_record_impacts()

    def _set_record_impacts(self):
        record_impacts = self._record_impacts
        if record_impacts is True:
            record_impacts = self.collimator_names
        if record_impacts is False or record_impacts is None:
            record_impacts = []
        if not hasattr(record_impacts, '__iter__') or isinstance(record_impacts, str):
            record_impacts = [record_impacts]
        if record_impacts and np.all([isinstance(self.line[elem], BaseCollimator) for elem in record_impacts]):
            for name in record_impacts:
                if name not in self.line.element_names:
                    raise ValueError(f"Trying to initialise impact table, but collimator {name} not found in line!")
            elements_to_record = list(set(record_impacts) - set(self._recording_elements))
            elements_to_stop = list(set(self._recording_elements) - set(record_impacts))
            if elements_to_stop:
                # These should have been stopped by the setter function above
                raise ValueError("Some elements are recording but are not supposed to!")
            # Any new elements that need recording but aren't recording yet?
            if elements_to_record:
                elements_to_record = [self.line[name] for name in elements_to_record]
                if self.impacts is None:
                    self._impacts = xt.start_internal_logging(io_buffer=self._io_buffer, capacity=self.capacity, \
                                                              elements=elements_to_record)
                else:
                    xt.start_internal_logging(io_buffer=self._io_buffer, capacity=self.capacity, \
                                              record=self.impacts, elements=elements_to_record)
                self._recording_elements = record_impacts
            self.impacts._coll_ids = {self.line.element_names.index(name): name for name in self._recording_elements}

    @property
    def capacity(self):
        return self._capacity

    @capacity.setter
    def capacity(self, capacity):
        capacity = int(capacity)
        if capacity < self.capacity:
            raise NotImplementedError("Shrinking of capacity not yet implemented!")
        elif capacity == self.capacity:
            return
        else:
            self._io_buffer.grow(capacity-self.capacity)
            if self.impacts is not None:
                # TODO: increase capacity of iobuffer AND of _impacts
                raise NotImplementedError

    @property
    def collimator_names(self):
        # TODO: should only sort whenever collimators are added to colldb
        # or when positions are updated
        db = self.colldb._colldb
#         names = list(db[db.active==True].index)
        names = list(db.index)
        if not np.any([s is None for s in db.s_center]):
            names.sort(key=lambda nn: db.loc[nn, 's_center'])
        return names

    def install_black_absorbers(self, names=None, *, verbose=False):
        line = self.line
        elements = []
        for name, gap, angle, length, side in zip(
                    self.collimator_names,
                    self.colldb.gap,
                    self.colldb.angle,
                    self.colldb.length,
                    self.colldb.side):
            if verbose: print(f"Installing {name:20} as BlackAbsorber")
            el = xc.BlackAbsorber(gap=gap, angle=angle, length=length, side=side, _tracking=False)

            # Check that collimator is not installed as different type
            # TODO: automatically replace collimator type and print warning
            if isinstance(line[name], tuple(collimator_classes)):
                raise ValueError(f"Trying to install {name} as {el.__class__.__name__},"
                               + f" but it is already installed as {line[name].__class__.__name__}!\n"
                               + f"Please reconstruct the line.")

            # TODO: only allow Marker elements, no Drifts!!
            #       How to do this with importing a line for MAD-X or SixTrack...?
            elif not isinstance(line[name], (xt.Marker, xt.Drift)) and not support_legacy_elements:
                raise ValueError(f"Trying to install {name} as {el.__class__.__name__},"
                               + f" but the line element to replace is not an xtrack.Marker "
                               + f"(or xtrack.Drift)!\nPlease check the name, or correct the "
                               + f"element.")
            elements.append(el)
        install_elements(line, self.collimator_names, elements, need_apertures=True)
        self._set_record_impacts()

    def install_everest_collimators(self, *, verbose=False):
        line = self.line
        elements = []
        for name, gap, angle, length, side, material, crystal, bending_radius, xdim, ydim, miscut, thick in zip(
                    self.collimator_names,
                    self.colldb.gap,
                    self.colldb.angle,
                    self.colldb.length,
                    self.colldb.side,
                    self.colldb.material,
                    self.colldb._colldb.crystal,
                    self.colldb._colldb.bending_radius,
                    self.colldb._colldb.xdim,
                    self.colldb._colldb.ydim,
                    self.colldb._colldb.miscut,
                    self.colldb._colldb.thick):
            mat = SixTrack_to_xcoll[material]
            if crystal is None:
                if verbose: print(f"Installing {name:20} as EverestCollimator")
                el = EverestCollimator(gap=gap, angle=angle, length=length, side=side, material=mat[0], _tracking=False)
            else:
                if verbose: print(f"Installing {name:20} as EverestCrystal")
                el = EverestCrystal(gap=gap, angle=angle, length=length, side=side, material=mat[1], lattice=crystal,
                                       bending_radius=bending_radius, xdim=xdim, ydim=ydim, miscut=miscut, thick=thick, _tracking=False)

            # Check that collimator is not installed as different type
            # TODO: automatically replace collimator type and print warning
            if isinstance(line[name], tuple(collimator_classes)):
                raise ValueError(f"Trying to install {name} as {el.__class__.__name__},"
                               + f" but it is already installed as {line[name].__class__.__name__}!\n"
                               + f"Please reconstruct the line.")

            # TODO: only allow Marker elements, no Drifts!!
            #       How to do this with importing a line for MAD-X or SixTrack...?
            elif not isinstance(line[name], (xt.Marker, xt.Drift)) and not support_legacy_elements:
                raise ValueError(f"Trying to install {name} as {el.__class__.__name__},"
                               + f" but the line element to replace is not an xtrack.Marker "
                               + f"(or xtrack.Drift)!\nPlease check the name, or correct the "
                               + f"element.")
            elements.append(el)
        install_elements(line, self.collimator_names, elements, need_apertures=True)
        self._set_record_impacts()

    def build_tracker(self, **kwargs):
        kwargs.setdefault('_buffer', self._buffer)
        kwargs.setdefault('io_buffer', self._io_buffer)
        if kwargs['_buffer'] != self._buffer:
            raise ValueError("Cannot build tracker with different buffer than the CollimatorManager buffer!")
        if kwargs['io_buffer'] != self._io_buffer:
            raise ValueError("Cannot build tracker with different io_buffer than the CollimatorManager io_buffer!")
        if '_context' in kwargs and kwargs['_context'] != self._buffer.context:
            raise ValueError("Cannot build tracker with different context than the CollimatorManager context!")
        self.line.build_tracker(**kwargs)
        self._set_record_impacts()


