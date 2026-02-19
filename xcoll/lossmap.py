# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import pandas as pd
from pathlib import Path
from types import GeneratorType
from concurrent.futures import ThreadPoolExecutor

import xtrack as xt
import xtrack.particles.pdg as pdg

from . import json
from .beam_elements import (collimator_classes, crystal_classes,
                            FlukaCollimator, FlukaCrystal,
                            Geant4Collimator, Geant4Crystal)
from .compare import deep_equal
from .general import __version__
from .plot import plot_lossmap, _resolve_zoom
from .constants import (USE_IN_LOSSMAP,
    LOST_ON_EVEREST_BLOCK, LOST_ON_EVEREST_COLL, LOST_ON_EVEREST_CRYSTAL,
    LOST_ON_FLUKA_BLOCK, LOST_ON_FLUKA_COLL, LOST_ON_FLUKA_CRYSTAL,
    LOST_ON_GEANT4_BLOCK, LOST_ON_GEANT4_COLL, LOST_ON_GEANT4_CRYSTAL,
    LOST_ON_BLACK_ABSORBER, LOST_ON_BLACK_CRYSTAL,)


class LossMap:
    _version_changes = ["0.6.0", "0.8.0.dev1+fluka"]

    def __init__(self, line=None, part=None, *, line_is_reversed=None, interpolation=None,
                 line_shift_s=0, weights=None, weight_function=None, verbose=True,
                 correct_aperture_absorption=True):
        self._line_is_reversed = None
        self._interpolation = None
        self._machine_length = None
        self._momentum = None
        self._aper_s = np.array([])
        self._aper_name = np.array([])
        self._aper_nabs = np.array([])
        self._aper_eabs = np.array([])
        self._aper_length = np.array([])
        self._aper_type = np.array([])
        self._aperbins = np.array([])
        self._aperbins_length = np.array([])
        self._aperbinned = np.array([])
        self._aperbinned_energy = np.array([])
        self._coll_s = np.array([])
        self._coll_name = np.array([])
        self._coll_nabs = np.array([])
        self._coll_eabs = np.array([])
        self._coll_length = np.array([])
        self._coll_type = np.array([])
        self._num_initial = 0
        self._tot_energy_initial = 0.
        self._beam_type = None
        self._cold_regions = None
        self._warm_regions = None
        self._s_range = {}
        self._xcoll = np.array([])
        self._date = np.array([])
        if part and line:
            self.add_particles(part=part, line=line, line_is_reversed=line_is_reversed,
                               interpolation=interpolation, line_shift_s=line_shift_s,
                               weights=weights, weight_function=weight_function,
                               correct_aperture_absorption=correct_aperture_absorption,
                               verbose=verbose)

    def __str__(self):
        return f"LossMap ({self.num_aperture_losses} losses on aperture and " \
             + f"{self.num_collimator_losses} losses on collimators - " \
             + f"{self.num_initial} initial {self.beam_type}s at {self.momentum:.3e} eV/c)"

    def __repr__(self):
        return f"<{str(self)} at {hex(id(self))}>"

    def __eq__(self, other):
        if not isinstance(other, LossMap):
            return False
        if not deep_equal(self.lossmap, other.lossmap):
            return False
        if not deep_equal(self.cold_regions, other.cold_regions):
            return False
        if not deep_equal(self.warm_regions, other.warm_regions):
            return False
        if not deep_equal(self.s_range, other.s_range):
            return False
        return True

    @classmethod
    def from_json(cls, file, verbose=True):
        lm = cls()
        lm.add_from_json(file, verbose=verbose)
        return lm

    def to_json(self, file):
        json.dump({
            'xcoll': self._xcoll,
            'momentum': self.momentum,
            'beam_type': self._beam_type,
            'num_initial': self.num_initial,
            'tot_energy_initial': self.tot_energy_initial,
            'cold_regions': self.cold_regions,
            'warm_regions': self.warm_regions,
            's_range': self.s_range,
            **self.lossmap,
            'date': self._date
        }, Path(file), indent=True)

    def save_summary(self, file):
        with open(Path(file), 'w') as fid:
            fid.write(self.summary.__repr__())


    # ================== #
    # === Properties === #
    # ================== #

    @property
    def num_losses(self):
        return self.num_aperture_losses + self.num_collimator_losses

    @property
    def num_aperture_losses(self):
        return int(self._aperbinned.sum() + self._aper_nabs.sum())

    @property
    def num_collimator_losses(self):
        return int(self._coll_nabs.sum())

    @property
    def num_initial(self):
        return self._num_initial

    @property
    def tot_energy(self):
        return self.tot_energy_aperture + self.tot_energy_collimator

    @property
    def tot_energy_collimator(self):
        return self._coll_eabs.sum()

    @property
    def tot_energy_aperture(self):
        return self._aperbinned_energy.sum()

    @property
    def tot_energy_initial(self):
        return self._tot_energy_initial


    @property
    def lossmap(self):
        coll_summary = self.summary[self.summary.n > 0].to_dict('list')
        coll_summary = {kk: np.array(vv) for kk, vv in coll_summary.items()}
        return {
                'machine_length':     self.machine_length,
                'interpolation':      self.interpolation,
                'reversed':           self.line_is_reversed,
                'collimator':         coll_summary,
                'aperture':           self.aperture_losses
            }

    @property
    def summary(self):
        return pd.DataFrame({
            'name':   self._coll_name,
            'n':      self._coll_nabs,
            'e':      self._coll_eabs,
            'length': self._coll_length,
            's':      self._coll_s,
            'type':   self._coll_type,
        }).sort_values("s")

    @property
    def aperture_losses(self):
        if self._interpolation:
            mask = self._aperbinned > 0
            aper_s = self._aperbins[:-1] + self.interpolation/2
            return {
                'idx_bins':    np.arange(len(self._aperbins)-1)[mask],
                's_bins':      aper_s[mask],
                'length_bins': self._aperbins_length[mask],
                'n_bins':      self._aperbinned[mask],
                'e_bins':      self._aperbinned_energy[mask]
            }
        else:
            return {
                'name':   self._aper_name,
                'n':      self._aper_nabs,
                'e':      self._aper_eabs,
                'length': self._aper_length,
                's':      self._aper_s,
                'type':   self._aper_type,
            }

    @property
    def interpolation(self):
        if self._interpolation == 0:
            return False
        return self._interpolation

    @interpolation.setter
    def interpolation(self, value):
        if value is None:
            if self._interpolation is not None:
                raise ValueError("Cannot unset interpolation once it has been set.")
        else:
            if value is False:
                value = 0
            elif value is True:
                value = 0.1
            elif value <= 0:
                raise ValueError("Interpolation step must be positive.")
            if self._interpolation is not None and not np.isclose(self._interpolation, value):
                raise ValueError("The interpolation step is different from the one "
                                "used to create the loss map.")
            if hasattr(value, 'tolist') and callable(value.tolist):
                value = value.tolist()
            self._interpolation = value

    @property
    def line_is_reversed(self):
        return self._line_is_reversed

    @line_is_reversed.setter
    def line_is_reversed(self, value):
        if value is not None:
            if self._line_is_reversed is not None:
                if self._line_is_reversed != value:
                    raise ValueError("The line_is_reversed is different from the one "
                                    "used to create the loss map.")
            if hasattr(value, 'tolist') and callable(value.tolist):
                value = value.tolist()
            self._line_is_reversed = value

    @property
    def momentum(self):
        return self._momentum

    @momentum.setter
    def momentum(self, value):
        if value is not None:
            if self._momentum is not None:
                if not np.isclose(self._momentum, value):
                    raise ValueError("The reference momentum is different from the one "
                                     "used to create the loss map.")
            if hasattr(value, 'tolist') and callable(value.tolist):
                value = value.tolist()
            self._momentum = value

    @property
    def beam_type(self):
        return pdg.get_name_from_pdg_id(self._beam_type)

    @beam_type.setter
    def beam_type(self, value):
        if value is not None:
            if self._beam_type is not None:
                if self._beam_type != value:
                    raise ValueError("The beam type is different from the one "
                                     "used to create the loss map.")
            if hasattr(value, 'tolist') and callable(value.tolist):
                value = value.tolist()
            if value == 0:
                value = 2212  # Assume proton if undefined
            self._beam_type = value

    @property
    def machine_length(self):
        return self._machine_length

    @machine_length.setter
    def machine_length(self, value):
        if value is not None:
            if self._machine_length is not None:
                if not np.isclose(self._machine_length, value):
                    raise ValueError("The machine_length is different from the one "
                                     "used to create the loss map.")
            if hasattr(value, 'tolist') and callable(value.tolist):
                value = value.tolist()
            self._machine_length = value

    @property
    def cold_regions(self):
        return self._cold_regions

    @cold_regions.setter
    def cold_regions(self, value):
        if value is not None:
            value = np.array(value)
            if self.cold_regions is not None:
                if self.warm_regions is not None:
                    raise ValueError("Cannot set both cold_regions and warm_regions.")
                if not deep_equal(self.cold_regions, value):
                    raise ValueError("New cold_regions do not match existing cold_regions.")
            self._cold_regions = value

    @property
    def warm_regions(self):
        return self._warm_regions

    @warm_regions.setter
    def warm_regions(self, value):
        if value is not None:
            value = np.array(value)
            if self.warm_regions is not None:
                if self.cold_regions is not None:
                    raise ValueError("Cannot set both cold_regions and warm_regions.")
                if not deep_equal(self.warm_regions, value):
                    raise ValueError("New warm_regions do not match existing warm_regions.")
            self._warm_regions = value

    @property
    def s_range(self):
        return self._s_range

    @s_range.setter
    def s_range(self, value):
        if value is not None:
            value = {kk: np.array(vv) for kk, vv in value.items()}
            if self._s_range:
                if set(value.keys()) != set(self._s_range.keys()):
                    raise ValueError("s_range is already set.")
                for kk, vv in value.items():
                    if not deep_equal(self._s_range[kk], vv):
                        raise ValueError("New s_range does not match existing s_range.")
            self._s_range = value


    # =============== #
    # === Methods === #
    # =============== #

    def update_metadata(self, metadata):
        if not isinstance(metadata, dict):
            raise ValueError("metadata must be a dictionary.")
        for kk, vv in metadata.items():
            if kk == 'cold_regions':
                self.cold_regions = vv
            elif kk == 'warm_regions':
                self.warm_regions = vv
            elif kk == 's_range':
                self.s_range = vv
            else:
                raise ValueError(f"Unknown metadata key '{kk}'.")


    def plot(self, **kwargs):
        cold_regions = self._cold_regions
        warm_regions = self._warm_regions
        xlim = kwargs.pop('xlim', None)
        if isinstance(xlim, str) and xlim in self.s_range:
            xlim = self.s_range[xlim]
        zoom = _resolve_zoom(kwargs.pop('zoom', None), allow_str=True)
        for zz in zoom:
            if isinstance(zz, str) and zz not in self.s_range:
                raise ValueError(f"Zoom string '{zz}' not found in `s_range`.")
        titles = kwargs.pop('titles', None)
        if zoom and np.all([isinstance(zz, str) for zz in zoom]):
            titles = ["Full ring",
                      *[f"Zoom on {zz} insertion" for zz in zoom]
                      ]
        zoom = [self.s_range[zz] if isinstance(zz, str) else zz for zz in zoom]
        energy = kwargs.pop('energy', None)
        if energy is None:
            energy = np.any([tt.startswith('Geant4') or tt.startswith('Fluka')
                             for tt in self._coll_type])
        beam_intensity = kwargs.pop('beam_intensity', None)
        if beam_intensity is not None and self.num_initial > 0:
            beam_intensity /= self.num_initial
        return plot_lossmap(self.lossmap, xlim=xlim, energy=energy, zoom=zoom,
                            cold_regions=cold_regions, warm_regions=warm_regions,
                            titles=titles, rescaled_beam_intensity=beam_intensity,
                            **kwargs)


    def add_particles(self, part, line, *, line_is_reversed=False, interpolation=None,
                      line_shift_s=0, weights=None, weight_function=None,
                      correct_aperture_absorption=True, verbose=True):
        """
        Add particles to the loss map. Aperture losses are interpolated and the
        collimator summary is updated.
        """
        # Check that collimators have been tracked
        geant4_coll = line.get_elements_of_type((Geant4Collimator, Geant4Crystal))[0]
        if len(geant4_coll) > 0 and np.all([coll._acc_ionisation_loss < 0 for coll in geant4_coll]):
            raise ValueError("Geant4Collimators have not been tracked, or LossMap already calculated")
        fluka_coll = line.get_elements_of_type((FlukaCollimator, FlukaCrystal))[0]
        if len(fluka_coll) > 0 and np.all([coll._acc_ionisation_loss < 0 for coll in fluka_coll]):
            raise ValueError("FlukaCollimators have not been tracked, or LossMap already calculated")
        if interpolation is not None:
            if self.interpolation is not None and not np.isclose(self.interpolation, interpolation):
                raise ValueError("The interpolation step is different from the one "
                                 "used to create the loss map.")
            self.interpolation = interpolation
        elif self.interpolation is None:
            self.interpolation = 0.1 # Default
        if not isinstance(correct_aperture_absorption, dict):
            coll_classes = list(set(collimator_classes) - set(crystal_classes))
            coll_elements = line.get_elements_of_type(coll_classes)[1]
            correct_aperture_absorption = {coll: correct_aperture_absorption
                                            for coll in coll_elements}
        for coll, this_aper_corr in correct_aperture_absorption.items():
            if this_aper_corr not in (True, False) \
            and this_aper_corr.lower() not in ('after', 'before', 'both'):
                raise ValueError("correct_aperture_absorption must be True, "
                                "False, 'after', 'before', 'both', or a dict "
                                "mapping collimator names to these values.")

        self.line_is_reversed = line_is_reversed
        self.machine_length = line.get_length()
        self.momentum = line.particle_ref.p0c[0]
        self.beam_type = line.particle_ref.pdg_id[0]
        if weights is None:
            if weight_function is None:
                weights = np.ones(len(part.x))
            else:
                weights = _create_weights_from_initial_state(part, weight_function)
        else:
            if weight_function is not None:
                raise ValueError("Use either 'weights' or 'weight_function', not both!")
            if len(weights) != len(part.x):
                raise ValueError("The length of the weights array must be equal to "
                                 "the number of particles.")

        if line_shift_s != 0:
            raise NotImplementedError("Line shift not implemented yet.")

        # Correct particles that are lost in aperture directly after collimator
        # (before interpolation to avoid issues with backtracking in collimators).
        aper_corr_coll = []
        for coll, this_aper_corr in correct_aperture_absorption.items():
            if this_aper_corr is True or isinstance(this_aper_corr, str) \
            and (this_aper_corr.lower() == 'both' or this_aper_corr.lower() == 'after'):
                aper_corr_coll.append(coll)
        self._correct_absorbed(part, line, verbose=verbose, aperture_loc='after', collimators=aper_corr_coll)

        # Loss location refinement
        if self._interpolation:
            self._interpolate(part, line, verbose=verbose)

        # Correct particles that are lost in aperture directly before collimator
        # (after interpolation to avoid moving too much losses incorrectly).
        aper_corr_coll = []
        for coll, this_aper_corr in correct_aperture_absorption.items():
            if (this_aper_corr is True and self._interpolation) or isinstance(this_aper_corr, str) \
            and (this_aper_corr.lower() == 'both' or this_aper_corr.lower() == 'before'):
                # Only correct before if interpolation was done (to avoid moving too much losses)
                aper_corr_coll.append(coll)
        self._correct_absorbed(part, line, verbose=verbose, aperture_loc='before', collimators=aper_corr_coll)

        self._make_coll_summary(part, line, line_shift_s, weights)
        self._get_aperture_losses(part, line, line_shift_s, weights)

        self._num_initial += np.sum((part.particle_id==part.parent_particle_id) & (part.particle_id>=0))
        self._tot_energy_initial += self.num_initial * self.momentum

        if 'collimation' in line.env.metadata:
            if 'cold_regions' in line.env.metadata['collimation']:
                self.cold_regions = line.env.metadata['collimation']['cold_regions']
            if 'warm_regions' in line.env.metadata['collimation']:
                self.warm_regions = line.env.metadata['collimation']['warm_regions']
            if 's_range' in line.env.metadata['collimation']:
                self.s_range = line.env.metadata['collimation']['s_range']
        if len(self._xcoll) == 0:
            self._xcoll = np.append(self._xcoll, __version__)
        elif __version__ not in self._xcoll:
            if verbose:
                print(f"Warning: xcoll version changed from {self._xcoll[-1]} to {__version__}.")
            self._xcoll = np.append(self._xcoll, __version__)
        self._date = np.append(self._date, pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'))


    def add_from_json(self, files, verbose=True):
        """
        Add loss map data from a JSON file or an iterable of files. The loss map
        is updated with the data from the file(s).
        """
        if not isinstance(files, (list, tuple, set, GeneratorType)):
            files = [files]
        files = list(files)

        # Define optimal parameters for parallelisation based on file size
        sz = Path(files[0]).stat().st_size
        if sz < 100_000:
            max_workers, chunksize = 16, 32
        elif sz < 500_000:
            max_workers, chunksize = 16, 16
        elif sz < 2_000_000:
            max_workers, chunksize = 16, 8
        else:
            max_workers, chunksize = 8, 4

        i = 0
        _xcoll = self._xcoll.tolist()
        _date = self._date.tolist()
        for lossmap in iter_lossmaps(files, max_workers=max_workers, chunksize=chunksize):
            LossMap._assert_valid_json(lossmap)
            if 'momentum' in lossmap:
                self.momentum = lossmap['momentum']
            if 'beam_type' in lossmap:
                self.beam_type = lossmap['beam_type']
            if 'xcoll' in lossmap:
                xcoll = [lossmap['xcoll']] if isinstance(lossmap['xcoll'], str) else lossmap['xcoll']
                _xcoll.extend(xcoll)
            if 'date' in lossmap:
                date = [lossmap['date']] if isinstance(lossmap['date'], str) else lossmap['date']
                _date.extend(date)
            self.line_is_reversed = lossmap['reversed']
            self.machine_length = lossmap['machine_length']
            self.interpolation = lossmap['interpolation']
            if 'cold_regions' in lossmap:
                self.cold_regions = lossmap['cold_regions']
            if 'warm_regions' in lossmap:
                self.warm_regions = lossmap['warm_regions']
            if 's_range' in lossmap:
                self.s_range = lossmap['s_range']
            self._load_coll_summary(lossmap['collimator'])
            self._load_aperture_losses(lossmap['aperture'])
            if 'num_initial' in lossmap and 'tot_energy_initial' in lossmap:
                self._num_initial += lossmap['num_initial']
                self._tot_energy_initial += lossmap['tot_energy_initial']
            elif not np.isclose(self.num_initial, 0) or not np.isclose(self.tot_energy_initial, 0):
                raise ValueError("num_initial and tot_energy_initial must be provided "
                                 "in the JSON file when adding to a non-empty LossMap.")
            i += 1
        if i == 0:
            raise ValueError("No valid files found.")
        _xcoll = np.unique(_xcoll)
        if verbose and len(_xcoll) > 1:
                print("Warning: Multiple xcoll versions are used in this loss map.")
        self._xcoll = np.unique(_xcoll)
        self._date = np.array(_date)
        if verbose:
            print(f"Loaded {i} file{'s' if i > 1 else ''} into loss map.")


    def _correct_absorbed(self, part, line, verbose, aperture_loc, collimators=[]):
        # Correct particles that are at an aperture directly before or after a collimator
        # TODO: should this be done if collimator has limited width/height?
        tt = line.get_table()
        for coll in collimators:
            elem = line.element_names.index(coll)
            if aperture_loc.lower() == 'before':
                aper = elem
                while not xt.line._is_aperture(line[aper], line):
                    aper = aper - 1 if aper > 0 else len(line.element_names) - 1
                mask  = (part.state == 0) & (part.at_element == aper)
                mask &= np.isclose(part.s, tt.rows[coll].s_start[0])
            elif aperture_loc.lower() == 'after':
                aper = elem
                while not xt.line._is_aperture(line[aper], line):
                    aper = aper + 1 if aper < len(line.element_names) - 1 else 0
                mask  = (part.state == 0) & (part.at_element == aper)
                mask &= np.isclose(part.s, tt.rows[coll].s_end[0])
            part.at_element[mask] = elem
            what_type = line[elem].__class__.__name__
            if what_type == 'EverestBlock':
                part.state[mask] = LOST_ON_EVEREST_BLOCK
            elif what_type == 'EverestCollimator':
                part.state[mask] = LOST_ON_EVEREST_COLL
            elif what_type == 'EverestCrystal':
                part.state[mask] = LOST_ON_EVEREST_CRYSTAL
            elif what_type == 'FlukaBlock':
                part.state[mask] = LOST_ON_FLUKA_BLOCK
            elif what_type == 'FlukaCollimator':
                part.state[mask] = LOST_ON_FLUKA_COLL
            elif what_type == 'FlukaCrystal':
                part.state[mask] = LOST_ON_FLUKA_CRYSTAL
            elif what_type == 'Geant4Block':
                part.state[mask] = LOST_ON_GEANT4_BLOCK
            elif what_type.startswith('Geant4Collimator'):
                part.state[mask] = LOST_ON_GEANT4_COLL
            elif what_type == 'Geant4Crystal':
                part.state[mask] = LOST_ON_GEANT4_CRYSTAL
            elif what_type == 'BlackAbsorber':
                part.state[mask] = LOST_ON_BLACK_ABSORBER
            elif what_type == 'BlackCrystal':
                part.state[mask] = LOST_ON_BLACK_CRYSTAL
            else:
                raise ValueError(f"Unknown collimator type {what_type}")
            if verbose and mask.sum() > 0:
                print(f"Found {mask.sum()} losses at {line.element_names[aper]}, "
                    + f"moved to {line.element_names[elem]}")


    def _interpolate(self, part, line, verbose=True):
        if len(part.s[part.state==0]) > 0:
            if verbose:
                print("Performing the aperture losses refinement.")
            time_dependent_was_enabled = False
            if line.enable_time_dependent_vars:
                line.enable_time_dependent_vars = False
                if verbose:
                    print("Temporarily disabled time-dependent variables in the line for loss location refinement.")
                time_dependent_was_enabled = True
            loss_loc_refinement = xt.LossLocationRefinement(
                line,
                n_theta = 360,            # Angular resolution
                r_max = 0.5,              # Maximum transverse aperture [m]
                dr = 50e-6,               # Transverse accuracy [m]
                ds = self.interpolation   # Longitudinal accuracy [m]
            )
            loss_loc_refinement.refine_loss_location(part)
            if time_dependent_was_enabled:
                line.enable_time_dependent_vars = True


    def _make_coll_summary(self, part, line, line_shift_s, weights):
        names = np.unique(line.get_elements_of_type(collimator_classes)[1])
        coll_mask = np.isin(part.state, USE_IN_LOSSMAP)
        coll_losses = np.array([line.element_names[i]
                                for i in part.at_element[coll_mask]])
        coll_lengths = [line[name].length for name in names]

        L = self.machine_length
        coll_pos = np.array([(line.get_s_position(name) + cl/2 + line_shift_s)%L
                    for name, cl in zip(names, coll_lengths)])
        if self.line_is_reversed:
            coll_pos = L - coll_pos

        coll_types = [line[name].__class__.__name__
                      for name in names]
        coll_weights = weights[coll_mask]
        nabs = [coll_weights[coll_losses == name].sum()
                for name in names]

        deposited_energy = {name: line[name]._acc_ionisation_loss
                                  if hasattr(line[name], '_acc_ionisation_loss')
                                  else 0.
                            for name in names}
        energy_weights = coll_weights * part.energy[coll_mask]
        eabs = [energy_weights[coll_losses == name].sum() + deposited_energy[name]
                for name in names]
        self._do_collimator_adding(coll_s=coll_pos, coll_name=names, coll_nabs=nabs,
                                   coll_eabs=eabs, coll_length=coll_lengths,
                                   coll_type=coll_types)

    def _load_coll_summary(self, colldata):
        coll_eabs = colldata['e'] if 'e' in colldata else np.zeros(len(colldata['s']))
        if 'type' in colldata:
            coll_types = colldata['type']
        else:
            coll_types = np.full(len(colldata['s']), "Unknown", dtype=object)
        self._do_collimator_adding(coll_s=colldata['s'], coll_name=colldata['name'],
                                   coll_nabs=colldata['n'], coll_eabs=coll_eabs,
                                   coll_length=colldata['length'], coll_type=coll_types)

    def _do_collimator_adding(self, coll_s, coll_name, coll_nabs, coll_eabs,
                              coll_length, coll_type):
        ds = 1e-6
        s_all      = np.concatenate([self._coll_s,    coll_s])
        nabs_all   = np.concatenate([self._coll_nabs, coll_nabs])
        eabs_all   = np.concatenate([self._coll_eabs, coll_eabs])
        name_all   = np.concatenate([self._coll_name,   coll_name])
        type_all   = np.concatenate([self._coll_type,   coll_type])
        length_all = np.concatenate([self._coll_length, coll_length])

        key = np.rint(s_all / ds).astype(np.int64)
        _, first, inv = np.unique(key, return_index=True, return_inverse=True)
        order = np.argsort(first)

        s_rep      = s_all[first]
        name_rep   = name_all[first]
        type_rep   = type_all[first]
        length_rep = length_all[first]

        _validate_str_meta(name_all, name_rep, inv, s_rep, "Collimator name", "coll_s")
        _validate_str_meta(type_all, type_rep, inv, s_rep, "Collimator type", "coll_s")
        _validate_float_meta(length_all, length_rep, inv, s_rep, "Collimator length", "coll_s")

        self._coll_s      = s_rep[order]
        self._coll_nabs   = np.bincount(inv, weights=nabs_all)[order]
        self._coll_eabs   = np.bincount(inv, weights=eabs_all)[order]
        self._coll_name   = name_rep[order]
        self._coll_type   = type_rep[order]
        self._coll_length = length_rep[order]


    def _make_aperture_bins(self):
        if self._aperbins.size == 0:
            L = self.machine_length
            end = L - np.mod(L, self.interpolation) + self.interpolation
            self._aperbins = np.linspace(0, end, int(np.ceil(L/self.interpolation)) + 1)
            # TODO: not entirely correct, as not all bins have the same lenght in practice: If
            # an active element like a magnet starts in the middle of the bin, that bin is
            # essentially shorter. Then the first aperture marker after the element represents
            # the particles that would be lost in this element, so there the aperture length is
            # equal to the element length. This is not true if the active element is a collimator
            # (then the losses are moved to the collimator anyway).
            # TODO: need length for non-interpolated case as well?
            self._aperbins_length = np.full(len(self._aperbins) - 1, self.interpolation, dtype=np.float64)
            self._aperbinned = np.zeros(len(self._aperbins) - 1, dtype=np.float64)
            self._aperbinned_energy = np.zeros(len(self._aperbins) - 1, dtype=np.float64)

    def _get_aperture_losses(self, part, line, line_shift_s, weights):
        aper_mask = part.state == 0
        if aper_mask.sum() == 0:
            return

        if self._interpolation:
            # Get s position per particle (lost on aperture)
            L = self.machine_length
            aper_s = np.mod(part.s[aper_mask] + line_shift_s, L)
            if self.line_is_reversed:
                aper_s = L - aper_s
            # Binned aperture losses
            self._make_aperture_bins()
            self._do_aperture_binning(aper_s=aper_s, aper_nabs=weights[aper_mask],
                                      aper_eabs=weights[aper_mask] * part.energy[aper_mask])
        else:
            aper_pos_dct, aper_length_dct = self._create_aperture_pos_dict(line, line_shift_s)
            part_idx, inv = np.unique(part.at_element[aper_mask], return_inverse=True)
            assert all([ii in aper_pos_dct for ii in part_idx])
            aper_nabs = np.bincount(inv, weights=weights[aper_mask])
            aper_eabs = np.bincount(inv, weights=weights[aper_mask] * part.energy[aper_mask])
            aper_pos = np.array([aper_pos_dct[ii] for ii in part_idx])
            aper_length = np.array([aper_length_dct[ii] for ii in part_idx])
            aper_name = np.array([line.element_names[idx] for idx in part_idx])
            aper_type = np.array([str(line.get(line.element_names[idx])) for idx in part_idx])
            self._do_aperture_adding(aper_name=aper_name,
                                     aper_nabs=aper_nabs,
                                     aper_eabs=aper_eabs,
                                     aper_length=aper_length,
                                     aper_s=aper_pos,
                                     aper_type=aper_type)

    def _load_aperture_losses(self, aperdata):
        if self._interpolation:
            self._make_aperture_bins()
            if not 'idx_bins' in aperdata and 's' in aperdata and 'n' in aperdata:
                # Old-style JSON file
                aper_eabs = aperdata['e'] if 'e' in aperdata else np.zeros(len(aperdata['s']))
                self._do_aperture_binning(aper_s=aperdata['s'], aper_nabs=aperdata['n'],
                                          aper_eabs=aper_eabs)
            else:
                self._aperbinned[aperdata['idx_bins']] += aperdata['n_bins']
                self._aperbinned_energy[aperdata['idx_bins']] += aperdata['e_bins']
        else:
            names = aperdata['name'] if 'name' in aperdata else np.array(['']*len(aperdata['s']))
            lengths = aperdata['length'] if 'length' in aperdata else np.ones(len(aperdata['s']))
            types = aperdata['type'] if 'type' in aperdata else np.array(['']*len(aperdata['s']))
            self._do_aperture_adding(aper_name=names,
                                     aper_nabs=aperdata['n'],
                                     aper_eabs=aperdata['e'],
                                     aper_length=lengths,
                                     aper_s=aperdata['s'],
                                     aper_type=types)

    def _do_aperture_binning(self, aper_s, aper_nabs, aper_eabs):
        binned = np.digitize(aper_s, bins=self._aperbins, right=False) - 1
        # We cannot directly use binned as an index, as repeated
        # indices would be ignored
        minlength = len(self._aperbins) - 1
        self._aperbinned += np.bincount(binned, weights=aper_nabs,
                                        minlength=minlength)
        self._aperbinned_energy += np.bincount(binned, weights=aper_eabs,
                                               minlength=minlength)

    def _create_aperture_pos_dict(self, line, line_shift_s):
        # Getting the position of aperture losses.
        # First, we identify the apertures before and after a given aperture, and
        # consider the mid-points between each of them and the current aperture.
        # The distance between these mid-points will act as the aperture's length,
        # while the position of the aperture will be given by the mid-point between
        # these two mid-points.
        # TODO: Need similar thing for interpolated case
        tt = line.get_table()
        mask = [ttt.startswith('Limit') for ttt in tt.element_type]
        mask_line = [ttt.startswith('Limit') for ttt in tt.element_type if ttt != '']
        assert len(mask_line) == len(line.element_names)
        tt = tt.rows[mask]
        idx = np.arange(len(mask_line))[mask_line]
        aper_prev = np.concatenate(([tt.s_start[-1]], tt.s_start[:-1]))
        aper_next = np.concatenate((tt.s_start[1:], [tt.s_start[0]]))
        aper_prev_mid = (aper_prev + tt.s_start) / 2
        aper_next_mid = (tt.s_start + aper_next) / 2
        L = self.machine_length
        pos = np.mod((aper_next_mid+aper_prev_mid)/2 + line_shift_s, L)
        if self._line_is_reversed:
            pos = L - pos
        pos = np.mod((aper_next_mid+aper_prev_mid)/2, line.get_length())
        aper_pos_dct    = dict(zip(idx, pos))
        aper_length_dct = dict(zip(idx, aper_next_mid-aper_prev_mid))
        return aper_pos_dct, aper_length_dct

    def _do_aperture_adding(self, aper_name, aper_nabs, aper_eabs, aper_length, aper_s, aper_type):
        ds = 1e-6
        s_all      = np.concatenate([self._aper_s,    aper_s])
        nabs_all   = np.concatenate([self._aper_nabs, aper_nabs])
        eabs_all   = np.concatenate([self._aper_eabs, aper_eabs])
        name_all   = np.concatenate([self._aper_name,   aper_name])
        type_all   = np.concatenate([self._aper_type,   aper_type])
        length_all = np.concatenate([self._aper_length, aper_length])

        key = np.rint(s_all / ds).astype(np.int64)
        _, first, inv = np.unique(key, return_index=True, return_inverse=True)
        order = np.argsort(first)

        s_rep      = s_all[first]
        name_rep   = name_all[first]
        type_rep   = type_all[first]
        length_rep = length_all[first]

        _validate_str_meta(name_all, name_rep, inv, s_rep, "Aperture name", "aper_s")
        _validate_str_meta(type_all, type_rep, inv, s_rep, "Aperture type", "aper_s")
        _validate_float_meta(length_all, length_rep, inv, s_rep, "Aperture length", "aper_s")

        self._aper_s      = s_rep[order]
        self._aper_nabs   = np.bincount(inv, weights=nabs_all)[order]
        self._aper_eabs   = np.bincount(inv, weights=eabs_all)[order]
        self._aper_name   = name_rep[order]
        self._aper_type   = type_rep[order]
        self._aper_length = length_rep[order]


    @staticmethod
    def _assert_valid_json(lossmap):
        # Get version number
        def get_ver(version):
            n_ver = sum([10**(3*(2-i))*int(j.split('rc')[0]) for i, j in enumerate(version.strip().split('.')[:3])])
            if 'rc' in version:
                n_ver += 0.5
            if len(version.strip().split('.')) > 3:
                n_ver += 0.2
            return n_ver
        if 'xcoll' in lossmap:
            vers = [get_ver(vv) for vv in lossmap['xcoll']]
            for vv in reversed(LossMap._version_changes):
                if any([v >= get_ver(vv) for v in vers]):
                    if not all([v >= get_ver(vv) for v in vers]):
                        raise ValueError("The xcoll versions in the JSON files are inconsistent.")
            n_ver = max(vers)
        else:
            n_ver = 0
        revision_1 = n_ver >= get_ver(LossMap._version_changes[0])
        revision_2 = n_ver >= get_ver(LossMap._version_changes[1])
        # General metadata
        if 'machine_length' not in lossmap:
            raise ValueError("The JSON file does not contain the machine length data.")
        if 'interpolation' not in lossmap:
            raise ValueError("The JSON file does not contain the interpolation data.")
        if 'reversed' not in lossmap:
            raise ValueError("The JSON file does not contain the reversed data.")
        if revision_1:
            if 'momentum' not in lossmap:
                raise ValueError("The JSON file does not contain the momentum data.")
        if revision_2:
            if 'beam_type' not in lossmap:
                raise ValueError("The JSON file does not contain the beam_type data.")
            if 'date' not in lossmap:
                raise ValueError("The JSON file does not contain the date data.")
            if 'cold_regions' not in lossmap:
                raise ValueError("The JSON file does not contain the cold_regions data.")
            if 'warm_regions' not in lossmap:
                raise ValueError("The JSON file does not contain the warm_regions data.")
            if 's_range' not in lossmap:
                raise ValueError("The JSON file does not contain the s_range data.")
            if 'num_initial' not in lossmap:
                raise ValueError("The JSON file does not contain the num_initial data.")
            if 'tot_energy_initial' not in lossmap:
                raise ValueError("The JSON file does not contain the tot_energy_initial data.")
        # Collimator losses
        if 'collimator' not in lossmap:
            raise ValueError("The JSON file does not contain the collimator data.")
        if 's' not in lossmap['collimator']:
            raise ValueError("The JSON file does not contain the collimator s data.")
        if 'name' not in lossmap['collimator']:
            raise ValueError("The JSON file does not contain the collimator name data.")
        if 'length' not in lossmap['collimator']:
            raise ValueError("The JSON file does not contain the collimator length data.")
        if 'n' not in lossmap['collimator']:
            raise ValueError("The JSON file does not contain the collimator n data.")
        if revision_1:
            if 'e' not in lossmap['collimator']:
                raise ValueError("The JSON file does not contain the collimator energy data.")
            if 'type' not in lossmap['collimator']:
                raise ValueError("The JSON file does not contain the collimator type data.")
        # Aperture losses
        if 'aperture' not in lossmap:
            raise ValueError("The JSON file does not contain the aperture data.")
        if revision_1 and lossmap['interpolation'] is not False:
            if 'idx_bins' not in lossmap['aperture']:
                raise ValueError("The JSON file does not contain the aperture idx_bins data.")
            if 'n_bins' not in lossmap['aperture']:
                raise ValueError("The JSON file does not contain the aperture n_bins data.")
            if 'e_bins' not in lossmap['aperture']:
                raise ValueError("The JSON file does not contain the aperture e_bins data.")
            if 'length_bins' not in lossmap['aperture']:
                raise ValueError("The JSON file does not contain the aperture length_bins data.")
        else:
            if 's' not in lossmap['aperture']:
                raise ValueError("The JSON file does not contain the aperture s data.")
            if 'n' not in lossmap['aperture']:
                raise ValueError("The JSON file does not contain the aperture n data.")
            if revision_1:
                if 'e' not in lossmap['aperture']:
                    raise ValueError("The JSON file does not contain the aperture e data.")
            if revision_2:
                if 'length' not in lossmap['aperture']:
                    raise ValueError("The JSON file does not contain the aperture length data.")
                if 'name' not in lossmap['aperture']:
                    raise ValueError("The JSON file does not contain the aperture name data.")
                if 'type' not in lossmap['aperture']:
                    raise ValueError("The JSON file does not contain the aperture type data.")


class MultiLossMap(LossMap):
    def __init__(self, *lms, **named_lms):
        """
        Create a MultiLossMap object from a list of LossMap objects.
        """
        super().__init__(None, None, line_is_reversed=None)
        del self._line_is_reversed
        self._lms = []
        self._lm_names = []
        if len(lms) > 0:
            if len(named_lms) > 0:
                raise ValueError("Specify the LossMaps as args or kwargs, not both.")
            if hasattr(lms[0], '__iter__'):
                if len(lms) > 1:
                    raise ValueError("Use args or a list of LossMap objects, not both.")
                if isinstance(lms[0], dict):
                    lms = lms[0].values()
                    self._lm_names = lms[0].keys()
                else:
                    lms = lms[0]
        elif len(named_lms) > 0:
            lms = named_lms.values()
            self._lm_names = named_lms.keys()
        for lm in lms:
            self.add_lossmap(lm)

    @property
    def lms(self):
        if len(self._lm_names) > 0:
            return dict(zip(self._lm_names, self._lms))
        else:
            return self._lms

    @property
    def line_is_reversed(self):
        raise AttributeError("The line_is_reversed property is not available "
                           + "for MultiLossMap objects.")

    @property
    def lossmap(self):
        coll_summary = self.summary[self.summary.n > 0].to_dict('list')
        coll_summary = {kk: np.array(vv) for kk, vv in coll_summary.items()}
        return {
                'collimator':      coll_summary,
                'aperture':        self.aperture_losses,
                'machine_length':  self.machine_length,
                'interpolation':   self.interpolation,
                'reversed':        self.line_is_reversed,
                'momentum':        self.momentum,
                'num_initial':     self.num_initial,
                'tot_energy_initial': self.tot_energy_initial
            }


    def add_lossmap(self, lm):
        """
        Add a LossMap object to the MultiLossMap object.
        """
        if not isinstance(lm, LossMap):
            raise ValueError("The input must be a LossMap object.")
        self._lms.append(lm)
        if self.machine_length is None:
            self._machine_length = lm.machine_length
        elif not np.isclose(self.machine_length, lm.machine_length):
            raise ValueError("The machine lengths of the loss maps are not the same.")
        if self.momentum is None:
            self._momentum = lm.momentum
        elif lm.momentum and not np.isclose(self.momentum, lm.momentum):
            raise ValueError("The momenta of the loss maps are not the same.")
        if self.interpolation is None:
            self._interpolation = lm.interpolation
        elif lm.interpolation is None and self.interpolation is not None:
            raise ValueError("The interpolations of the loss maps are not the same.")
        elif not np.isclose(self.interpolation, lm.interpolation):
            raise ValueError("The interpolations of the loss maps are not the same.")
        if self.cold_regions is None:
            self._cold_regions = lm.cold_regions
        elif not deep_equal(self.cold_regions, lm.cold_regions):
            raise ValueError("The cold_regions of the loss maps are not the same.")
        if self.warm_regions is None:
            self._warm_regions = lm.warm_regions
        elif not deep_equal(self.warm_regions, lm.warm_regions):
            raise ValueError("The warm_regions of the loss maps are not the same.")
        if self.s_range is None:
            self._s_range = lm.s_range
        elif not deep_equal(self.s_range, lm.s_range):
            raise ValueError("The s_range of the loss maps are not the same.")
        self._num_initial += lm.num_initial
        self._tot_energy_initial += lm.tot_energy_initial
        self._do_collimator_adding(coll_s=lm._coll_s, coll_name=lm._coll_name,
                                   coll_nabs=lm._coll_nabs, coll_eabs=lm._coll_eabs,
                                   coll_length=lm._coll_length, coll_type=lm._coll_type)
        if self.interpolation:
            if len(lm._aperbins) == 0:
                # No aperture in this loss map
                return
            elif len(self._aperbins) == 0:
                self._aperbins = lm._aperbins
                self._aperbins_length = lm._aperbins_length
                self._aperbinned = lm._aperbinned
                self._aperbinned_energy = lm._aperbinned_energy
            else:
                if len(self._aperbins) != len(lm._aperbins) \
                or not deep_equal(self._aperbins, lm._aperbins):
                    raise ValueError("The number of bins of the loss maps are not the same.")
                self._aperbinned += lm._aperbinned
                self._aperbinned_energy += lm._aperbinned_energy
        else:
            self._do_aperture_adding(aper_s=lm._aper_s, aper_nabs=lm._aper_nabs,
                                     aper_eabs=lm._aper_eabs, aper_length=lm._aper_length,
                                     aper_name=lm._aper_name, aper_type=lm._aper_type)


def _create_weights_from_initial_state(part, function):
    if len(function) == 4:
        return function[0](part.x)*function[1](part.px)*\
               function[2](part.y)*function[3](part.py)
    elif len(function) == 6:
        return function[0](part.x)*function[1](part.px)*\
               function[2](part.y)*function[3](part.py)*\
               function[4](part.zeta)*function[5](part.delta)
    else:
        raise NotImplementedError


def _validate_str_meta(values_all, values_rep, inv, s_rep, label, label_s):
    bad = values_all != values_rep[inv]
    if not np.any(bad):
        return
    i = np.flatnonzero(bad)[0]
    g = inv[i]
    vals = np.unique(values_all[inv == g])
    raise ValueError(
        f"{label} mismatch at {label_s} ≈ {s_rep[g]!r}:\n"
        f"  representative: {values_rep[g]!r}\n"
        f"  values found: {vals}"
    )

def _validate_float_meta(values_all, values_rep, inv, s_rep, label,
                         label_s, atol=1e-12, rtol=1e-12):
    ok = np.isclose(values_all, values_rep[inv], atol=atol, rtol=rtol) | (
        np.isnan(values_all) & np.isnan(values_rep[inv])
    )
    if np.all(ok):
        return
    i = np.flatnonzero(~ok)[0]
    g = inv[i]
    vals = np.unique(values_all[inv == g])
    raise ValueError(
        f"{label} mismatch at {label_s} ≈ {s_rep[g]!r}:\n"
        f"  representative: {values_rep[g]!r}\n"
        f"  values found: {vals}"
    )


def iter_lossmaps(paths, max_workers=16, chunksize=16):
    # Parallel loading of JSON files (validation and aggregation is not parallellised to avoid race conditions)
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        yield from ex.map(json.load, paths, chunksize=chunksize)
