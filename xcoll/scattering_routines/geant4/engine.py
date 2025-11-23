# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025                  #
# ######################################### #

import os
import sys
import numpy as np
from pathlib import Path
from numbers import Number

from .rpyc import launch_rpyc_with_port # remove after geant4 bugfix

import xobjects as xo

from .std_redirect import pin_python_stdio
from .bdsim_config import create_bdsim_config_file, get_collimators_from_input_file
from ..engine import BaseEngine
from ...general import _pkg_root
try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath


class Geant4Engine(BaseEngine):

    _xofields = BaseEngine._xofields | {
        '_relative_energy_cut':        xo.Float64
        # 'random_freeze_state':         xo.Int64,  # to be implemented; number of randoms already sampled, such that this can be taken up again later
    }

    _int32 = True
    _uses_input_file = True
    _uses_run_folder = True

    _depends_on = [BaseEngine]

    def __init__(self, **kwargs):
        # Set element classes dynamically
        from ...beam_elements import Geant4Collimator, Geant4Crystal
        self.__class__._element_classes = (Geant4Collimator, Geant4Crystal)
        # Initialise geant4-only defaults
        self._g4link = None
        self._server = None # remove after geant4 bugfix
        self._conn = None # remove after geant4 bugfix
        super().__init__(**kwargs)
        self._already_started = False
        # Set default values for new properties
        self.relative_energy_cut = None
        self.reentry_protection_enabled = None


    # ======================
    # === New Properties ===
    # ======================

    @property
    def relative_energy_cut(self):
        return self._relative_energy_cut

    @relative_energy_cut.setter
    def relative_energy_cut(self, val):
        if val is None:
            val = 0.1
        if not isinstance(val, Number) or val <= 0:
            raise ValueError("`relative_energy_cut` has to be a strictly postive number!")
        self._relative_energy_cut = val

    @property
    def reentry_protection_enabled(self):
        return self._reentry_protection_enabled

    @reentry_protection_enabled.setter
    def reentry_protection_enabled(self, val):
        if val is False:
            print("Warning: Disabling re-entry protection can lead to crashes!")
            print("         Only disable if you know what you are doing.")
        elif val is None:
            try:
                import rpyc
            except ImportError as e:
                val = False
            else:
                val = True
        elif not isinstance(val, bool):
            raise ValueError("`reentry_protection_enabled` has to be a boolean!")
        self._reentry_protection_enabled = val

    @property
    def lower_momentum_cut(self):
        return self._lower_momentum_cut

    @lower_momentum_cut.setter
    def lower_momentum_cut(self, val):
        if val is None:
            val = 0.0  # TODO: keep in mind that relative_energy_cut must be consistent with this
        else:
            raise NotImplementedError  # TODO
        if not isinstance(val, Number) or val < 0:
            raise ValueError("`lower_momentum_cut` has to be a non-negative number!")
        self._lower_momentum_cut = val

    @property
    def photon_lower_momentum_cut(self):
        return self._photon_lower_momentum_cut

    @photon_lower_momentum_cut.setter
    def photon_lower_momentum_cut(self, val):
        if val is None:
            val = 0.0  # TODO: keep in mind that relative_energy_cut must be consistent with this
        else:
            raise NotImplementedError  # TODO
        if not isinstance(val, Number) or val < 0:
            raise ValueError("`photon_lower_momentum_cut` has to be a non-negative number!")
        self._photon_lower_momentum_cut = val

    @property
    def electron_lower_momentum_cut(self):
        return self._electron_lower_momentum_cut

    @electron_lower_momentum_cut.setter
    def electron_lower_momentum_cut(self, val):
        if val is None:
            val = 0.0  # TODO: keep in mind that relative_energy_cut must be consistent with this
        else:
            raise NotImplementedError  # TODO
        if not isinstance(val, Number) or val < 0:
            raise ValueError("`electron_lower_momentum_cut` has to be a non-negative number!")
        self._electron_lower_momentum_cut = val


    # ============================
    # === Overwrite Properties ===
    # ============================

    @property
    def capacity(self):
        return None  # Geant4 capacity is dynamic

    @property
    def relative_capacity(self):
        return None  # Geant4 capacity is dynamic


    # =================================
    # === Base methods to overwrite ===
    # =================================

    def _set_engine_properties(self, **kwargs):
        self._set_property('relative_energy_cut', kwargs)
        kwargs = super()._set_engine_properties(**kwargs)
        return kwargs

    def _pre_input(self, **kwargs):
        coll_id = 1
        for el in self._element_dict.values():
            el.geant4_id = f'XcollG4.{coll_id}'  # TODO: will be provided by new BDSIM interface
            coll_id += 1
        return kwargs

    def _generate_input_file(self, **kwargs):
        input_file, kwargs = create_bdsim_config_file(element_dict=self._element_dict,
                                particle_ref=self.particle_ref, verbose=self.verbose,
                                **kwargs)
        # The only thing left in kwargs are parameters to start the engine
        return input_file, kwargs

    def _start_engine(self, **kwargs):
        from ...beam_elements import BaseCrystal

        Ekin = self.particle_ref.energy0 - self.particle_ref.mass0
        try:
            from g4interface import XtrackInterface
        except (ModuleNotFoundError, ImportError) as error:
            self.stop(clean=True)
            self._warn(error)
            return

        if self.reentry_protection_enabled:
            ### remove this part after geant4 bug fixed
            try:
                import rpyc
            except (ModuleNotFoundError, ImportError) as e:
                raise ImportError("Failed to import rpyc. Cannot connect to BDSIM.") from e
            self._server, port = launch_rpyc_with_port()
            self._conn = rpyc.classic.connect('localhost', port=port)
            self._conn._config['sync_request_timeout'] = 1240 # Set timeout to 1240 seconds
            self._conn.execute('import sys')
            self._conn.execute(f'sys.path.append("{self._environment.data_dir.as_posix()}")')  # For g4interface on server
            self._conn.execute(f'sys.path.append("{(_pkg_root / "scattering_routines" / "geant4").as_posix()}")')
            self._conn.execute('import engine_server')
            self._g4link = self._conn.namespace['engine_server'].BDSIMServer()
            self._g4link.XtrackInterface(bdsimConfigFile='geant4_input.gmad',
                                         referencePdgId=self.particle_ref.pdg_id,
                                         referenceEk=Ekin,
                                         relativeEnergyCut=self.relative_energy_cut,
                                         seed=self.seed, batchMode=True)
        else:
            if self._already_started:
                self.stop(clean=True)
                raise RuntimeError("Cannot restart Geant4 engine in non-reentry-safe mode. "
                                 + "Please exit this Python process. Do pip install rpyc "
                                 + "to avoid this limitation.")

            # Take iostream copies before we construct XtrackInterface (i.e. before FDRedirect runs)
            # to avoid python output being redirected by FDRedirect in C
            with pin_python_stdio():
                self._g4link = XtrackInterface(bdsimConfigFile='geant4_input.gmad',
                                               referencePdgId=self.particle_ref.pdg_id,
                                               referenceEk=Ekin,
                                               relativeEnergyCut=self.relative_energy_cut,
                                               seed=self.seed, batchMode=True)

        for el in self._element_dict.values():
            side = 2 if el._side == -1 else el._side
            jaw_L = 0.1 if el.jaw_L is None else el.jaw_L
            jaw_R = -0.1 if el.jaw_R is None else el.jaw_R
            tilt_L = 0.0 if el.tilt_L is None else el.tilt_L
            tilt_R = 0.0 if el.tilt_R is None else el.tilt_R
            jaw_L -= 1.e-9  # Correct for 1e-9 shift that is added in BDSIM
            jaw_R += 1.e-9  # Correct for 1e-9 shift that is added in BDSIM
            if jaw_L <= 0:
                self.stop(clean=True)
                raise NotImplementedError(f"Geant4Collimator {el.name} has negative left jaw: "
                               + f"jaw_L={el.jaw_L}! BDSIM cannot handle this.")
            if jaw_R >= 0:
                self.stop(clean=True)
                raise NotImplementedError(f"Geant4Collimator {el.name} has positive right jaw: "
                               + f"jaw_R={el.jaw_R}! BDSIM cannot handle this.")
            self._g4link.addCollimator(f'{el.geant4_id}', el.material.geant4_name, el.length,
                                    apertureLeft=jaw_L, apertureRight=-jaw_R,
                                    rotation=np.deg2rad(el.angle),
                                    xOffset=0, yOffset=0, side=side,
                                    jawTiltLeft=tilt_L, jawTiltRight=tilt_R,
                                    isACrystal=isinstance(el, BaseCrystal))
        self._already_started = True

    def _stop_engine(self, **kwargs):
        self._g4link = None
        if self.reentry_protection_enabled and self._server: # remove after geant4 bugfix
            self._server.terminate() # remove after geant4 bugfix
            self._server = None # remove after geant4 bugfix
        return kwargs

    def _is_running(self):
        return self._g4link is not None

    def _get_input_files_to_clean(self, input_file, cwd, **kwargs):
        if cwd is None or input_file is None:
            return []
        return [cwd / input_file]

    def _get_output_files_to_clean(self, input_file, cwd, **kwargs):
        if cwd is None:
            return []
        files_to_delete = ['rpyc.log', 'geant4.out', 'geant4.err',
                           'engine.out', 'engine.err', 'root.out',
                           'root.err']
        return [cwd / f for f in files_to_delete]

    def _match_input_file(self):
        # Read the elements in the input file and compare to the elements in the engine,
        # overwriting parameters where necessary
        input_dict = get_collimators_from_input_file(self.input_file)
        for name in input_dict:
            if name not in self._element_dict:
                raise ValueError(f"Element {name} in input file not found in engine!")
        for name, ee in self._element_dict.items():
            if name not in input_dict:
                self._print(f"Warning: Geant4Collimator {name} not in Geant4 input file! "
                          + f"Maybe it was fully open. Deactivated")
                self._deactivate_element(ee)
                continue
            self._assert_element(ee)
            ee.geant4_id = input_dict[name]['geant4_id']
            if not np.isclose(ee.length, input_dict[name]['length'], atol=1e-9):
                self._print(f"Warning: Length of {name} differs from input file "
                        + f"({ee.length} vs {input_dict[name]['length']})! Overwritten.")
                ee.length = input_dict[name]['length']
            if not np.isclose(ee.angle, input_dict[name]['angle'], atol=1e-9):
                self._print(f"Warning: Angle of {name} differs from input file "
                        + f"({ee.angle} vs {input_dict[name]['angle']})! Overwritten.")
                ee.angle = input_dict[name]['angle']
            if ee.material.geant4_name != input_dict[name]['material'] \
            and ee.material.name != input_dict[name]['material']:
                raise ValueError(f"Material of {name} differs from input file "
                            + f"({ee.material.geant4_name or ee.material.name} "
                            + f"vs {input_dict[name]['material']})!")
            jaw = input_dict[name]['jaw']
            if jaw is not None and not hasattr(jaw, '__iter__'):
                jaw = [jaw, -jaw]
            if jaw is None or (jaw[0] is None and jaw[1] is None):
                ee.jaw = None
            else:
                if jaw[0] is None:
                    if ee.side != 'right':
                        self._print(f"Warning: {name} is right-sided in the input file, "
                                + "but not in the line! Overwritten by the former.")
                        ee.side = 'right'
                elif ee.jaw_L is None or not np.isclose(ee.jaw_L, jaw[0], atol=1e-9):
                    self._print(f"Warning: Jaw_L of {name} differs from input file "
                            + f"({ee.jaw_L} vs {jaw[0]})! Overwritten.")
                    ee.jaw_L = jaw[0]
                if jaw[1] is None:
                    if ee.side != 'left':
                        self._print(f"Warning: {name} is left-sided in the input file, "
                                + f"but not in the line! Overwritten by the former.")
                        ee.side = 'left'
                elif ee.jaw_R is None or not np.isclose(ee.jaw_R, jaw[1], atol=1e-9):
                    self._print(f"Warning: Jaw_R of {name} differs from input file "
                            + f"({ee.jaw_R} vs {jaw[1]})! Overwritten.")
                    ee.jaw_R = jaw[1]
            tilt = input_dict[name]['tilt']
            if not hasattr(tilt, '__iter__'):
                tilt = [tilt, -tilt]
            if ee.side != 'right' and not np.isclose(ee.tilt_L, tilt[0], atol=1e-9):
                self._print(f"Warning: Tilt_L of {name} differs from input file "
                        + f"({ee.tilt_L} vs {tilt[0]})! Overwritten.")
                ee.tilt_L = tilt[0]
            if ee.side != 'left' and not np.isclose(ee.tilt_R, tilt[1], atol=1e-9):
                self._print(f"Warning: Tilt_R of {name} differs from input file "
                        + f"({ee.tilt_R} vs {tilt[1]})! Overwritten.")
                ee.tilt_R = tilt[1]
