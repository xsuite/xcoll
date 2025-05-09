# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import pandas as pd
from pathlib import Path
import json

import xtrack as xt
import xobjects as xo

from .beam_elements import collimator_classes, crystal_classes
from .general import __version__
from .plot import _plot_lossmap_base


class LossMap:
    def __init__(self, line, part, *, line_is_reversed, interpolation=None,
                 line_shift_s=0, weights=None, weight_function=None, verbose=True):
        self._line_is_reversed = line_is_reversed
        self._machine_length = line.get_length() if line else None
        self._interpolation = interpolation
        self._momentum = line.particle_ref.p0c[0] if line else None
        self._aper_s = np.array([])
        self._aper_nabs = np.array([])
        self._aper_eabs = np.array([])
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
        self._xcoll = __version__
        self._energy_waring_given = False
        if part and line:
            self.add_particles(part=part, line=line, line_shift_s=line_shift_s,
                               weights=weights, weight_function=weight_function,
                               verbose=verbose)

    def __str__(self):
        return f"LossMap ({self.num_aperture_losses} losses on aperture and " \
             + f"{self.num_collimator_losses} losses on collimators)"

    def __repr__(self):
        return f"<{str(self)} at {hex(id(self))}>"

    @classmethod
    def from_json(cls, file, verbose=True):
        lm = cls(None, None, line_is_reversed=None)
        lm.add_from_json(file, verbose=verbose)
        return lm

    def to_json(self, file):
        with open(Path(file), 'w') as fid:
            json.dump({
                'xcoll': self._xcoll,
                'date':  pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
                **self.lossmap
            }, fid, indent=True, cls=xo.JEncoder)

    def save_summary(self, file):
        with open(Path(file), 'w') as fid:
            fid.write(self.summary.__repr__())


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
                'momentum':        self._momentum,
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
                'idx_bins': np.arange(len(self._aperbins)-1)[mask],
                's_bins': aper_s[mask],
                'n_bins': self._aperbinned[mask],
                'e_bins': self._aperbinned_energy[mask],
                'length_bins': self._aperbins_length[mask],
            }
        else:
            return {
                's': self._aper_s,
                'n': self._aper_nabs,
                'e': self._aper_eabs
            }

    @property
    def line_is_reversed(self):
        return self._line_is_reversed

    @property
    def machine_length(self):
        return self._machine_length

    @property
    def interpolation(self):
        return self._interpolation

    @property
    def momentum(self):
        return self._momentum

    @property
    def num_losses(self):
        return self.num_aperture_losses + self.num_collimator_losses

    @property
    def num_aperture_losses(self):
        return int(self._aperbinned.sum() + self._aper_nabs.sum())

    @property
    def num_collimator_losses(self):
        return int(self._coll_nabs.sum())


    def plot(self, *, norm="total", ax=None, xlim=None, ylim=None, legend=True,
             grid=True, energy=False, show=True, savefig=None):
        return _plot_lossmap_base(self.lossmap, norm=norm, ax=ax, xlim=xlim,
                            ylim=ylim, legend=legend, grid=grid, energy=energy,
                            show=show, savefig=savefig)

    def add_particles(self, part, line, *, line_shift_s=0, weights=None,
                      weight_function=None, verbose=True):
        """
        Add particles to the loss map. Aperture losses are interpolated and the
        collimator summary is updated.
        """
        if self.interpolation is None:
            self._interpolation = 0.1
        if self._machine_length:
            if not np.isclose(self._machine_length, line.get_length()):
                raise ValueError("The line length is different from the one used "
                                 "to create the loss map.")
        else:
            self._machine_length = line.get_length()
        if self._momentum:
            if not np.isclose(self._momentum, line.particle_ref.p0c[0]):
                raise ValueError("The reference momentum is different from the one "
                                 "used to create the loss map.")
        else:
            self._momentum = line.particle_ref.p0c[0]
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

        # Correct particles that are lost in aperture directly after collimator.
        # These should be absorbed.
        self._correct_absorbed(part, line, verbose=verbose)

        # Loss location refinement
        if self._interpolation:
            self._interpolate(part, line, verbose=verbose)

        self._make_coll_summary(part, line, line_shift_s, weights)
        self._get_aperture_losses(part, line_shift_s, weights)


    def add_from_json(self, file, verbose=True):
        """
        Add loss map data from a JSON file or an iterable of files. The loss map
        is updated with the data from the file(s).
        """
        if hasattr(file, '__iter__') and not isinstance(file, (str, bytes, Path)):
            i = 0
            for ff in file:
                self.add_from_json(ff)
                i += 1
            if verbose:
                print(f"Loaded {i} files into loss map.")
            return
        with open(Path(file), 'r') as fid:
            lossmap = json.load(fid)
        new_style_json = LossMap._assert_valid_json(lossmap)
        if new_style_json:
            if self._momentum is None:
                self._momentum = lossmap['momentum']
            elif not np.isclose(self._momentum, lossmap['momentum']):
                raise ValueError("The momentum is different from the one used "
                                 "to create the loss map.")
            if self._xcoll != lossmap['xcoll'] and verbose:
                print("Warning: The xcoll version is different from the one used "
                    "to create the loss map.")
        if self.line_is_reversed is None:
            self._line_is_reversed = lossmap['reversed']
        elif self.line_is_reversed != lossmap['reversed']:
            raise ValueError("The line_is_reversed is different from the one used "
                             "to create the loss map.")
        if self.machine_length is None:
            self._machine_length = lossmap['machine_length']
        elif not np.isclose(self.machine_length, lossmap['machine_length']):
            raise ValueError("The line length is different from the one used "
                             "to create the loss map.")
        if self.interpolation is None:
            self._interpolation = lossmap['interpolation']
        elif not np.isclose(self.interpolation, lossmap['interpolation']):
            raise ValueError("The interpolation step is different from the one used "
                             "to create the loss map.")
        self._load_coll_summary(lossmap['collimator'])
        self._load_aperture_losses(lossmap['aperture'])


    def _correct_absorbed(self, part, line, verbose=True):
        # Correct particles that are at an aperture directly after a collimator
        coll_classes = list(set(collimator_classes) - set(crystal_classes))
        coll_elements = line.get_elements_of_type(coll_classes)[1]
        for idx, elem in enumerate(part.at_element):
            if part.state[idx] == 0:
                if elem == 0:
                    prev_elem = len(line.element_names) - 1
                else:
                    prev_elem = elem - 1
                if line.element_names[prev_elem] in coll_elements:
                    if verbose:
                        print(f"Found at {line.element_names[elem]}, "
                            + f"moved to {line.element_names[elem-1]}")
                    part.at_element[idx] = elem - 1
                    what_type = line[elem-1].__class__.__name__
                    if what_type == 'EverestCollimator':
                        part.state[idx] = -331
                    elif what_type == 'EverestCrystal':
                        part.state[idx] = -332
                    elif what_type == 'FlukaCollimator':
                        part.state[idx] = -334   # TODO: what if crystal?
                    elif what_type == 'Geant4Collimator':
                        part.state[idx] = -337   # TODO: what if crystal?
                    elif what_type == 'BlackAbsorber':
                        part.state[idx] = -340
                    else:
                        raise ValueError(f"Unknown collimator type {what_type}")


    def _interpolate(self, part, line, verbose=True):
        if len(part.s[part.state==0]) > 0:
            if verbose:
                print("Performing the aperture losses refinement.")
            loss_loc_refinement = xt.LossLocationRefinement(
                line,
                n_theta = 360,            # Angular resolution
                r_max = 0.5,              # Maximum transverse aperture [m]
                dr = 50e-6,               # Transverse accuracy [m]
                ds = self.interpolation   # Longitudinal accuracy [m]
            )
            loss_loc_refinement.refine_loss_location(part)


    def _make_coll_summary(self, part, line, line_shift_s, weights):
        names = line.get_elements_of_type(collimator_classes)[1]
        coll_mask = (part.state <= -330) & (part.state >= -350)
        coll_losses = np.array([line.element_names[i]
                                  for i in part.at_element[coll_mask]])
        coll_lengths = [line[name].length for name in names]

        L = self.machine_length
        coll_pos = [(line.get_s_position(name) + cl/2 + line_shift_s)%L
                    for name, cl in zip(names, coll_lengths)]
        if self.line_is_reversed:
            coll_pos = L - coll_pos

        coll_types = [line[name].__class__.__name__
                      for name in names]
        coll_weights = weights[coll_mask]
        nabs = [coll_weights[coll_losses == name].sum()
                for name in names]
        energy_weights = coll_weights * part.energy[coll_mask]
        eabs = [energy_weights[coll_losses == name].sum()
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
        # TODO: this can be done smarter, masks instead of loops, though
        # one should keep in mind that indices can be repeated and then
        # masking does not work correctly
        for ss, nn, nabs, eabs, ll, tt in zip(coll_s, coll_name, coll_nabs,
                                        coll_eabs, coll_length, coll_type):
            idx, = np.where(self._coll_name == nn)
            if len(idx) == 0:
                self._coll_s = np.append(self._coll_s, ss)
                self._coll_name = np.append(self._coll_name, nn)
                self._coll_nabs = np.append(self._coll_nabs, nabs)
                self._coll_eabs = np.append(self._coll_eabs, eabs)
                self._coll_length = np.append(self._coll_length, ll)
                self._coll_type = np.append(self._coll_type, tt)
            else:
                idx = idx[0]
                if not np.isclose(self._coll_length[idx], ll):
                    raise ValueError(f"Length of {nn} is different from the one used "
                                    + "to create the loss map.")
                if self._coll_type[idx] != tt and tt != 'Unknown':
                    raise ValueError(f"Type of {nn} is different from the one used to "
                                    + "create the loss map.")
                if not np.isclose(self._coll_s[idx], ss):
                    raise ValueError(f"s position of {nn} is different from the one "
                                    + "used to create the loss map.")
                self._coll_nabs[idx] += nabs
                self._coll_eabs[idx] += eabs


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

    def _get_aperture_losses(self, part, line_shift_s, weights):
        aper_mask = part.state == 0
        if len(part.s[aper_mask]) == 0:
            return

        # Get s position per particle (lost on aperture)
        L = self.machine_length
        aper_s = np.mod(part.s[aper_mask] + line_shift_s, L)
        if self._line_is_reversed:
            aper_s = L - aper_s

        if self._interpolation:
            # Binned aperture losses
            self._make_aperture_bins()
            self._do_aperture_binning(aper_s=aper_s, aper_nabs=weights[aper_mask],
                                      aper_eabs=weights[aper_mask] * part.energy[aper_mask])
        else:
            # Aperture losses at exact s positions (because no interpolation performed)
            # TODO: need correct lengths to scale
            aper_pos       = np.unique(aper_s)
            aper_weights   = weights[aper_mask]
            aper_nabs      = [aper_weights[aper_s == j].sum() for j in aper_pos]
            energy_weights = aper_weights * part.energy[aper_mask]
            aper_eabs      = [energy_weights[aper_s == s].sum() for s in aper_pos]
            self._do_aperture_adding(aper_s=aper_pos, aper_nabs=aper_nabs,
                                     aper_eabs=aper_eabs)

    def _load_aperture_losses(self, aperdata):
        if self._interpolation:
            self._make_aperture_bins()
            if not 'idx_bins' in aperdata and 's' in aperdata and 'n' in aperdata:
                # Old-style JSON file
                aper_eabs = aperdata['e'] if 'e' in aperdata else np.zeros(len(aperdata['s']))
                self._do_aperture_binning(aper_s=aperdata['s'], aper_nabs=aperdata['n'],
                                          aper_eabs=aper_eabs)
            else:
                self._aperbinned[aperdata['idx_bins']] = aperdata['n_bins']
                self._aperbinned_energy[aperdata['idx_bins']] = aperdata['e_bins']
        else:
            self._do_aperture_adding(aper_s=aperdata['s'], aper_nabs=aperdata['n'],
                                     aper_eabs=aperdata['e'])

    def _do_aperture_binning(self, aper_s, aper_nabs, aper_eabs):
        binned = np.digitize(aper_s, bins=self._aperbins, right=False) - 1
        # We cannot directly use binned as an index, as repeated
        # indices would be ignored
        minlength = len(self._aperbins) - 1
        self._aperbinned += np.bincount(binned, weights=aper_nabs,
                                        minlength=minlength)
        self._aperbinned_energy += np.bincount(binned, weights=aper_eabs,
                                               minlength=minlength)

    def _do_aperture_adding(self, aper_s, aper_nabs, aper_eabs):
        # TODO: this can be done smarter, masks instead of loops, though
        # one should keep in mind that indices can be repeated and then
        # masking does not work correctly
        for ss, nabs, eabs, ll in zip(aper_s, aper_nabs, aper_eabs):
            idx, = np.where(self._aper_s == ss)
            if len(idx) == 0:
                self._aper_s = np.append(self._aper_s, ss)
                self._aper_nabs = np.append(self._aper_nabs, nabs)
                self._aper_eabs = np.append(self._aper_eabs, eabs)
            else:
                idx = idx[0]
                self._aper_nabs[idx] += nabs
                self._aper_eabs[idx] += eabs


    @staticmethod
    def _assert_valid_json(lossmap):
        is_new_format = 'xcoll' in lossmap
        if 'machine_length' not in lossmap:
            raise ValueError("The JSON file does not contain the machine length data.")
        if 'interpolation' not in lossmap:
            raise ValueError("The JSON file does not contain the interpolation data.")
        if 'reversed' not in lossmap:
            raise ValueError("The JSON file does not contain the reversed data.")
        if is_new_format:
            if 'momentum' not in lossmap:
                raise ValueError("The JSON file does not contain the momentum data.")
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
        if is_new_format:
            if 'e' not in lossmap['collimator']:
                raise ValueError("The JSON file does not contain the collimator energy data.")
            if 'type' not in lossmap['collimator']:
                raise ValueError("The JSON file does not contain the collimator type data.")
        if 'aperture' not in lossmap:
            raise ValueError("The JSON file does not contain the aperture data.")
        if is_new_format and lossmap['interpolation'] is not None:
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
            if is_new_format:
                if 'e' not in lossmap['aperture']:
                    raise ValueError("The JSON file does not contain the aperture e data.")
                if 'length' not in lossmap['aperture']:
                    raise ValueError("The JSON file does not contain the aperture length data.")
        return is_new_format


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
                'momentum':        self._momentum,
            }


    def add_lossmap(self, lm):
        """
        Add a LossMap object to the MultiLossMap object.
        """
        if not isinstance(lm, LossMap):
            raise ValueError("The input must be a LossMap object.")
        self._lms.append(lm)
        if self._machine_length is None:
            self._machine_length = lm.machine_length
        elif not np.isclose(self._machine_length, lm.machine_length):
            raise ValueError("The machine lengths of the loss maps are not the same.")
        if self._momentum is None:
            self._momentum = lm.momentum
        elif lm.momentum and not np.isclose(self._momentum, lm.momentum):
            raise ValueError("The momenta of the loss maps are not the same.")
        if self._interpolation is None:
            self._interpolation = lm.interpolation
        elif not np.isclose(self._interpolation, lm.interpolation):
            raise ValueError("The interpolations of the loss maps are not the same.")
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
                or not np.allclose(self._aperbins, lm._aperbins):
                    raise ValueError("The number of bins of the loss maps are not the same.")
                self._aperbinned += lm._aperbinned
                self._aperbinned_energy += lm._aperbinned_energy
        else:
            self._do_aperture_adding(aper_s=lm._aper_s, aper_nabs=lm._aper_nabs,
                                     aper_eabs=lm._aper_eabs)


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

