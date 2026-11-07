# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025                  #
# ######################################### #

import gc
import numpy as np
from numbers import Number

import xobjects as xo

from .rpyc import launch_rpyc_with_port # remove after geant4 bugfix
from .bdsim_config import create_bdsim_config_file, get_collimators_from_input_file
from .reference_masses import geant4_masses_src
from ..engine import BaseEngine
from ...general import _pkg_root


class Geant4Engine(BaseEngine):

    _xofields = BaseEngine._xofields #| {
    #     random_freeze_state':         xo.Int64,  # to be implemented; number of randoms already sampled, such that this can be taken up again later
    # }

    _int32 = True
    _uses_input_file = True
    _uses_run_folder = True
    # These cannot (yet) be set on the Geant4 interface
    _physics_settings_veto_list = [
        'hadron_lower_momentum_cut',
        'photon_lower_momentum_cut',
        'electron_lower_momentum_cut',
        'include_showers',
        'include_single_coulomb',
        'include_multiple_coulomb',
        'include_elastic',
        'include_inelastic',
        'include_ionisation_fluctuations',
        'include_pair_production',
        'include_bremsstrahlung'
    ]

    _depends_on = [BaseEngine]

    _extra_c_sources = [geant4_masses_src]

    def __init__(self, **kwargs):
        # Set element classes dynamically
        from ...beam_elements import Geant4Collimator, Geant4CollimatorTip, Geant4Crystal
        self.__class__._element_classes = (Geant4Collimator, Geant4CollimatorTip, Geant4Crystal)
        # Initialise geant4-only defaults
        self._g4link = None
        self._server = None # remove after geant4 bugfix
        self._conn = None # remove after geant4 bugfix
        super().__init__(**kwargs)
        self._already_started = False
        # Set default values for new properties
        self.reentry_protection_enabled = None


    # ======================
    # === New Properties ===
    # ======================

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


    # =================================
    # === Base methods to overwrite ===
    # =================================

    def _set_engine_properties(self, **kwargs):
        kwargs = super()._set_engine_properties(**kwargs)
        self._set_property('reentry_protection_enabled', kwargs)
        return kwargs


    def _reset_engine_settings(self):
        self.reentry_protection_enabled = None


    def _pre_input(self, **kwargs):
        coll_id = 1
        for el in self._element_dict.values():
            el.geant4_id = f'XcollG4.{coll_id}'  # TODO: will be provided by new BDSIM interface
            coll_id += 1
        return kwargs


    def _generate_input_file(self, **kwargs):
        input_file, kwargs = create_bdsim_config_file(element_dict=self._element_dict,
                                particle_ref=self.particle_ref, verbose=self.verbose,
                                cwd=self.cwd, **kwargs)
        # The only thing left in kwargs are parameters to start the engine
        return input_file, kwargs


    def _start_engine(self, **kwargs):
        from ...beam_elements import BaseCrystal, Geant4CollimatorTip

        Ekin = self.particle_ref.energy0[0] - self.particle_ref.mass0
        try:
            from g4interface import XtrackInterface
        except (ModuleNotFoundError, ImportError) as error:
            self.stop()
            self._warn(error)
            return

        # An interface built against a BDSIM without the tip taper (bdsim-collaboration/bdsim#228), or
        # an older build of it, would silently drop these settings: a flat jaw in a
        # 2 m wide box instead of the requested geometry
        if not getattr(XtrackInterface, 'supportsTipTaper', lambda: False)():
            for el in self._element_dict.values():
                if isinstance(el, Geant4CollimatorTip) and any(
                        getattr(el, field) != 0
                        for field in ('taper_angle', 'jaw_depth', 'taper_depth', 'jaw_width')):
                    self.stop(clean=True)
                    raise RuntimeError(
                        f"Geant4CollimatorTip {el.name} uses taper_angle, jaw_depth, taper_depth "
                        f"or jaw_width, which the compiled Xcoll-BDSIM interface does not support! "
                        f"This needs a BDSIM that includes bdsim-collaboration/bdsim#228 "
                        f"(https://github.com/bdsim-collaboration/bdsim/pull/228), and the "
                        f"interface to be compiled against it: xc.geant4.interface.compile().")

        if self.reentry_protection_enabled:
            ### remove this part after geant4 bug fixed
            try:
                import rpyc
            except (ModuleNotFoundError, ImportError) as e:
                self.stop()
                raise ImportError("Failed to import rpyc. Cannot connect to BDSIM.") from e
            self._server, port = launch_rpyc_with_port(log_path=self.cwd / "rpyc.log")
            self._conn = rpyc.classic.connect('localhost', port=port)
            self._conn._config['sync_request_timeout'] = 1240 # Set timeout to 1240 seconds
            self._conn.execute('import sys')
            self._conn.execute(f'sys.path.append("{self.interface.lib_dir.as_posix()}")')  # For g4interface on server
            self._conn.execute(f'sys.path.append("{(_pkg_root / "scattering_routines" / "geant4").as_posix()}")')
            self._conn.execute('import engine_server')
            self._g4link = self._conn.namespace['engine_server'].BDSIMServer()
            self._g4link.XtrackInterface(bdsimConfigFile=self.input_file.as_posix(),
                                         referencePdgId=self.particle_ref.pdg_id[0],
                                         referenceEk=Ekin,
                                         relativeEnergyCut=self._physics_settings.relative_energy_cut,
                                         seed=self.seed, batchMode=True,
                                         workdir=self.cwd.as_posix())
        else:
            if self._already_started:
                self.stop(clean=True)
                raise RuntimeError("Cannot restart Geant4 engine in non-reentry-safe mode. "
                                 + "Please exit this Python process. Do pip install rpyc "
                                 + "to avoid this limitation.")

            self._g4link = XtrackInterface(bdsimConfigFile=self.input_file.as_posix(),
                                           referencePdgId=self.particle_ref.pdg_id[0],
                                           referenceEk=Ekin,
                                           relativeEnergyCut=self._physics_settings.relative_energy_cut,
                                           seed=self.seed, batchMode=True,
                                           workdir=self.cwd.as_posix())

        for el in self._element_dict.values():
            side = 2 if el._side == -1 else el._side
            jaw_L = 0.1 if el.jaw_L is None else el.jaw_L
            jaw_R = -0.1 if el.jaw_R is None else el.jaw_R
            tilt_L = 0.0 if el.tilt_L is None else el.tilt_L
            tilt_R = 0.0 if el.tilt_R is None else el.tilt_R
            jaw_L -= 1.e-9  # Correct for 1e-9 shift that is added in BDSIM
            jaw_R += 1.e-9  # Correct for 1e-9 shift that is added in BDSIM
            if jaw_L < jaw_R:
                self.stop(clean=True)
                raise ValueError(f"Geant4Collimator {el.name} has jaw_L < jaw_R: "
                             + f"jaw_L={el.jaw_L}, jaw_R={el.jaw_R}!")
            xOffset = 0
            yOffset = 0
            if jaw_L <= 0:
                xOffset = (jaw_L + jaw_R)/2
                jaw_L -= xOffset
                jaw_R -= xOffset
                xOffset_temp = xOffset*np.cos(np.deg2rad(el.angle_L)) - yOffset*np.sin(np.deg2rad(el.angle_L))
                yOffset = xOffset*np.sin(np.deg2rad(el.angle_L)) + yOffset*np.cos(np.deg2rad(el.angle_L))
                xOffset = xOffset_temp
            if jaw_R >= 0:
                xOffset = (jaw_L + jaw_R)/2
                jaw_L -= xOffset
                jaw_R -= xOffset
                xOffset_temp = xOffset*np.cos(np.deg2rad(el.angle_R)) - yOffset*np.sin(np.deg2rad(el.angle_R))
                yOffset = xOffset*np.sin(np.deg2rad(el.angle_R)) + yOffset*np.cos(np.deg2rad(el.angle_R))
                xOffset = xOffset_temp
            tip_material = ''
            tip_thickness = 0
            taper_angle = 0
            taper_depth = 0
            jaw_half_height = 0
            horizontal_width = 0
            if isinstance(el, Geant4CollimatorTip):
                tip_material = el.tip_material.geant4_name
                tip_thickness = el.tip_thickness
                taper_angle = el.taper_angle
                taper_depth = el.taper_depth
                jaw_half_height = el.jaw_width/2
                if el.jaw_depth > 0 and side == 0 \
                and not np.isclose(jaw_L, -jaw_R, rtol=0, atol=1e-12):
                    # Centre the jaws on the offset, so that a single horizontalWidth
                    # gives both of them exactly jaw_depth
                    shift = (jaw_L + jaw_R)/2
                    jaw_L -= shift
                    jaw_R -= shift
                    xOffset += shift*np.cos(np.deg2rad(el.angle))
                    yOffset += shift*np.sin(np.deg2rad(el.angle))
                try:
                    horizontal_width = self._tip_horizontal_width(el, jaw_L, jaw_R, side,
                                                                  tilt_L, tilt_R)
                except ValueError:
                    self.stop(clean=True)
                    raise
            self._g4link.addCollimator(f'{el.geant4_id}', el.material.geant4_name,
                                       tip_material, tip_thickness, el.length,
                                       apertureLeft=jaw_L, apertureRight=-jaw_R,
                                       rotation=np.deg2rad(el.angle),
                                       xOffset=xOffset, yOffset=yOffset, side=side,
                                       jawTiltLeft=tilt_L, jawTiltRight=tilt_R,
                                       taperAngle=taper_angle,
                                       horizontalWidth=horizontal_width,
                                       taperDepth=taper_depth,
                                       jawHalfHeight=jaw_half_height,
                                       isACrystal=isinstance(el, BaseCrystal))
        self._already_started = True

    @staticmethod
    def _tip_horizontal_width(el, jaw_L, jaw_R, side, tilt_L, tilt_R):
        # Returns the horizontalWidth to give BDSIM (0 means its default of 2 m).
        # As in BDSCollimatorJaw, each jaw spans from its edge out to
        # horizontalWidth/2. With a taper the jaw is `length` long at the beam
        # and grows with depth up to taper_depth, if set, or else over the
        # whole jaw (see taper_extension).
        # Whatever BDSIM would reject is checked here, as its link API does not
        # raise but silently replaces a failing collimator by a drift.
        half_gaps = []
        if side in (0, 1):
            half_gaps.append(jaw_L)
        if side in (0, 2):
            half_gaps.append(-jaw_R)
        if el.jaw_depth < 0:
            raise ValueError(f"Geant4CollimatorTip {el.name}: jaw_depth cannot be negative!")
        if el.jaw_width < 0:
            raise ValueError(f"Geant4CollimatorTip {el.name}: jaw_width cannot be negative!")
        if el.jaw_depth > 0:
            if min(half_gaps) < 0:
                raise ValueError(f"Geant4CollimatorTip {el.name}: jaw_depth is not supported for "
                               + "a jaw that crosses the beam axis!")
            horizontal_width = 2*(max(half_gaps) + el.jaw_depth)
        else:
            horizontal_width = 2.0
        depth = min(horizontal_width/2 - gap for gap in half_gaps)
        if depth <= el.tip_thickness:
            raise ValueError(f"Geant4CollimatorTip {el.name}: jaw depth ({depth*1e3:.3f} mm) "
                           + f"must be larger than tip_thickness ({el.tip_thickness*1e3:.3f} mm)!")
        if el.taper_angle != 0:
            if not 0 < el.taper_angle < np.pi/2:
                raise ValueError(f"Geant4CollimatorTip {el.name}: taper_angle must be in "
                               + "(0, 90) degrees!")
            if tilt_L != 0 or tilt_R != 0:
                raise ValueError(f"Geant4CollimatorTip {el.name}: a taper cannot be combined "
                               + "with a jaw tilt!")
            if el.taper_depth < 0:
                raise ValueError(f"Geant4CollimatorTip {el.name}: taper_depth cannot be negative!")
            if el.jaw_depth <= 0 and el.taper_depth <= 0:
                raise ValueError(f"Geant4CollimatorTip {el.name}: a taper needs jaw_depth or "
                               + "taper_depth, otherwise the wedge would span BDSIM's 2 m wide "
                               + "box and make the jaws several metres long!")
        return horizontal_width if el.jaw_depth > 0 else 0

    def _stop_engine(self, **kwargs):
        del self._g4link
        gc.collect()
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
            from ...beam_elements import Geant4CollimatorTip
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
            if isinstance(ee, Geant4CollimatorTip) or 'tip_material' in input_dict[name]:
                if not isinstance(ee, Geant4CollimatorTip):
                    raise ValueError(f"Element {name} is not a Geant4CollimatorTip "
                                    + "in the line, but it has tip material in the input file!")
                if 'tip_material' not in input_dict[name] or 'tip_thickness' not in input_dict[name]:
                    raise ValueError(f"Element {name} is a Geant4CollimatorTip, "
                                    + "but it has no tip material in the input file!")
                if ee.tip_material.geant4_name != input_dict[name]['tip_material']:
                    self._print(f"Warning: Tip material of {name} differs from input file "
                            + f"({ee.tip_material.geant4_name} vs {input_dict[name]['tip_material']})! Overwritten.")
                    ee.tip_material.geant4_name = input_dict[name]['tip_material']
                if not np.isclose(ee.length, input_dict[name]['tip_thickness'], atol=1e-9):
                    self._print(f"Warning: Tip thickness of {name} differs from input file "
                            + f"({ee.tip_thickness} vs {input_dict[name]['tip_thickness']})! Overwritten.")
                    ee.tip_thickness = input_dict[name]['tip_thickness']
                taper_angle = input_dict[name].get('taper_angle', 0)
                if not np.isclose(ee.taper_angle, taper_angle, atol=1e-9):
                    self._print(f"Warning: Taper angle of {name} differs from input file "
                            + f"({ee.taper_angle} vs {taper_angle})! Overwritten.")
                    ee.taper_angle = taper_angle
                jaw_depth = input_dict[name].get('jaw_depth', 0)
                if not np.isclose(ee.jaw_depth, jaw_depth, atol=1e-9):
                    self._print(f"Warning: Jaw depth of {name} differs from input file "
                            + f"({ee.jaw_depth} vs {jaw_depth})! Overwritten.")
                    ee.jaw_depth = jaw_depth
                taper_depth = input_dict[name].get('taper_depth', 0)
                if not np.isclose(ee.taper_depth, taper_depth, atol=1e-9):
                    self._print(f"Warning: Taper depth of {name} differs from input file "
                            + f"({ee.taper_depth} vs {taper_depth})! Overwritten.")
                    ee.taper_depth = taper_depth
                jaw_width = input_dict[name].get('jaw_width', 0)
                if not np.isclose(ee.jaw_width, jaw_width, atol=1e-9):
                    self._print(f"Warning: Jaw width of {name} differs from input file "
                            + f"({ee.jaw_width} vs {jaw_width})! Overwritten.")
                    ee.jaw_width = jaw_width
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
