# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025                  #
# ######################################### #

import numpy as np
from pathlib import Path
from numbers import Number
from contextlib import redirect_stdout, redirect_stderr

from subprocess import Popen # remove after geant4 bugfix
import socket # remove after geant4 bugfix
import time # remove after geant4 bugfix
from .rpyc import launch_rpyc_with_port # remove after geant4 bugfix

import xobjects as xo

from ..engine import BaseEngine
from ...general import _pkg_root
try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath


class Geant4Engine(BaseEngine):

    _xofields = {**BaseEngine._xofields,
        '_relative_energy_cut': xo.Float64,
        '_bdsim_config_file':   xo.String
#         'random_freeze_state':    xo.Int64,  # to be implemented; number of randoms already sampled, such that this can be taken up again later
    }

    _int32 = True
    # _uses_input_file = True
    _uses_run_folder = True

    def __init__(self, **kwargs):
        # Set element classes dynamically
        from ...beam_elements import Geant4Collimator, Geant4Crystal
        self.__class__._element_classes = (Geant4Collimator, Geant4Crystal)
        # Initialise geant4-only defaults
        self._g4link = None
        self._server = None # remove after geant4 bugfix
        self._conn = None # remove after geant4 bugfix
        kwargs['_bdsim_config_file'] = ''.ljust(256)
        super().__init__(**kwargs)
        self.relative_energy_cut = None # To set default value
        self.bdsim_config_file = None   # To set default value


    # ======================
    # === New Properties ===
    # ======================

    @property
    def relative_energy_cut(self):
        return self._relative_energy_cut

    @relative_energy_cut.setter
    def relative_energy_cut(self, val):
        if val is None:
            val = 0.15
        if not isinstance(val, Number) or val <= 0:
            raise ValueError("`relative_energy_cut` has to be a strictly postive number!")
        self._relative_energy_cut = val

    @property
    def bdsim_config_file(self):
        if self._bdsim_config_file == '':
            return None
        return FsPath(self._bdsim_config_file)

    @bdsim_config_file.setter
    def bdsim_config_file(self, val):
        if val is None:
            self._bdsim_config_file = ''
        else:
            if not isinstance(val, (str,Path)):
                raise ValueError("`bdsim_config_file` has to be a string or Path!")
            self._bdsim_config_file = FsPath(val).expanduser().resolve().as_posix()


    # =================================
    # === Base methods to overwrite ===
    # =================================

    def _set_engine_properties(self, **kwargs):
        kwargs = super()._set_engine_properties(**kwargs)
        self._set_property('relative_energy_cut', kwargs)
        self._set_property('bdsim_config_file', kwargs)
        return kwargs

    def _use_input_file(self, input_file=None, **kwargs):
        # Temporary hack to set geant4 IDs at the correct moment (before the setting of attributes is locked,
        # but after the engine has assigned the element_dict). When the gmad input file is auto-generated,
        # that script can also set the IDs.
        coll_id = 1
        for el in self._element_dict.values():
            el.geant4_id = coll_id # TODO: will be provided by new BDSIM interface
            coll_id += 1
        return kwargs

    def _start_engine(self, **kwargs):
        from ...beam_elements import BaseCrystal
        if self.bdsim_config_file is None:
            raise ValueError("`bdsim_config_file` must be set before starting the Geant4 engine!")

        Ekin = self.particle_ref.energy0 - self.particle_ref.mass0
        pdg_id = self.particle_ref.pdg_id

        ### revert after geant4 bug fixed
        ### try:
        ###     import collimasim as cs
        ### except ImportError as e:
        ###     raise ImportError("Failed to import collimasim. Cannot connect to BDSIM.")
        ### else:
        ###     self._g4link = cs.XtrackInterface(bdsimConfigFile=self.bdsim_config_file.as_posix(),
        ###                                       referencePdgId=self.particle_ref.pdg_id,
        ###                                       referenceEk=Ekin / 1e9, ### BDSIM expects GeV
        ###                                       relativeEnergyCut=self.relative_energy_cut,
        ###                                       seed=self.seed, batchMode=True)

        ### remove the following lines after geant4 bug fixed
        import rpyc
        self._server, port = launch_rpyc_with_port(log_path="rpyc.out")
        self._conn = rpyc.classic.connect('localhost', port=port)
        self._conn._config['sync_request_timeout'] = 1240 # Set timeout to 1240 seconds
        self._conn.execute('import sys')
        self._conn.execute(f'sys.path.append("{(_pkg_root / "scattering_routines" / "geant4").as_posix()}")')
        self._conn.execute('import engine_server')
        self._conn.execute('import collimasim as cs')
        self._g4link = self._conn.namespace['engine_server'].BDSIMServer()
        self._g4link.XtrackInterface(bdsimConfigFile=self.bdsim_config_file.as_posix(),
                                    referencePdgId=self.particle_ref.pdg_id,
                                    referenceEk=Ekin / 1e9, # BDSIM expects GeV
                                    relativeEnergyCut=self.relative_energy_cut,
                                    seed=self.seed, batchMode=True)
        ### remove down to here after geant4 bug fixed

        for el in self._element_dict.values():
            side = 2 if el._side == -1 else el._side
            jaw_L = 0.1 if el.jaw_L is None else el.jaw_L
            jaw_R = -0.1 if el.jaw_R is None else el.jaw_R
            tilt_L = 0.0 if el.tilt_L is None else el.tilt_L
            tilt_R = 0.0 if el.tilt_R is None else el.tilt_R
            # TODO: should geant4_id be a string or an int?
            self._g4link.addCollimator(f'{el.geant4_id}', el.material, el.length,
                                    apertureLeft=jaw_L-1.e-9,  # Correct for 1e-9 shift that is added in BDSIM
                                    apertureRight=-jaw_R-1.e-9,
                                    rotation=np.deg2rad(el.angle),
                                    xOffset=0, yOffset=0, side=side,
                                    jawTiltLeft=tilt_L, jawTiltRight=tilt_R,
                                    isACrystal=isinstance(el, BaseCrystal))

    def _stop_engine(self, **kwargs):
        self._g4link = None
        if self._server: # remove after geant4 bugfix
            self._server.terminate() # remove after geant4 bugfix
            self._server = None # remove after geant4 bugfix
        return kwargs

    def _is_running(self):
        return self._g4link is not None
