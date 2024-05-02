# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import json
from pathlib import Path
import numpy as np
import pandas as pd

from .beam_elements import BaseCollimator, BlackAbsorber, EverestCollimator, EverestCrystal, collimator_classes, element_classes
from .install import install_elements
from .colldb import CollimatorDatabase
from .interaction_record import InteractionRecord
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

    _init_vars = ['_colldb', 'line', 'beam']
    _init_var_defaults = {'_colldb': None, 'beam': None}

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

    def __getitem__(self, name):
        return self.colldb[name]

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
            el.emittance = self.colldb.emittance
            elements.append(el)
        install_elements(line, self.collimator_names, elements, need_apertures=True)

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
            el.emittance = self.colldb.emittance
            elements.append(el)
        install_elements(line, self.collimator_names, elements, need_apertures=True)


