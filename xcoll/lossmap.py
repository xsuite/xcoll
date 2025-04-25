# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pandas as pd
from pathlib import Path
import json

import xtrack as xt
import xobjects as xo

from .beam_elements import collimator_classes, crystal_classes
from .general import __version__


class LossMap:
    def __init__(self, line, part, *, line_is_reversed, interpolation=None,
                 line_shift_s=0, weights=None, weight_function=None, verbose=True):
        self._line_is_reversed = line_is_reversed
        self._machine_length = line.get_length() if line else None
        self._interpolation = interpolation
        self._momentum = line.particle_ref.p0c[0] if line else None
        self._aperture = {'s': [], 'name': [], 'n': [], 'e': []}
        self._summary = {'s': [], 'name': [], 'n': [], 'e': [], 'length': [],
                         'type': []}
        self._xcoll = __version__
        self._energy_waring_given = False
        if part and line:
            self.add_particles(part=part, line=line, line_shift_s=line_shift_s,
                               weights=weights, weight_function=weight_function,
                               verbose=verbose)

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
            fid.write(self._summary.__repr__())


    @property
    def lossmap(self):
        coll_summary = self.summary[self.summary.n > 0].to_dict('list')
        return {
                'collimator':     coll_summary,
                'aperture':       self._aperture,
                'machine_length': self.machine_length,
                'interpolation':  self.interpolation,
                'reversed':       self.line_is_reversed,
                'momentum':       self._momentum,
            }

    @property
    def summary(self):
        return pd.DataFrame({
            'name':   self._summary['name'],
            'n':      self._summary['n'],
            'e':      self._summary['e'],
            'length': self._summary['length'],
            's':      self._summary['s'],
            'type':   self._summary['type'],
        })

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

        self._add_summary(self._make_coll_summary(part, line, line_shift_s, weights))
        self._add_aperture(self._get_aperture_losses(part, line, line_shift_s, weights),
                           verbose=verbose)


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
        LossMap._assert_valid_json(lossmap)
        if 'momentum' in lossmap:
            if self._momentum is None:
                self._momentum = lossmap['momentum']
            elif not np.isclose(self._momentum, lossmap['momentum']):
                raise ValueError("The momentum is different from the one used "
                                 "to create the loss map.")
        if 'xcoll' in lossmap and self._xcoll != lossmap['xcoll'] and verbose:
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
        self._add_summary(lossmap['collimator'])
        self._add_aperture(lossmap['aperture'])

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
        collimator_names = line.get_elements_of_type(collimator_classes)[1]
        coll_mask = (part.state <= -330) & (part.state >= -350)
        coll_losses = np.array([line.element_names[i]
                                  for i in part.at_element[coll_mask]])
        coll_lengths = [line[name].length for name in collimator_names]

        L = self.machine_length
        coll_pos = [(line.get_s_position(name) + cl/2 + line_shift_s)%L
                    for name, cl in zip(collimator_names, coll_lengths)]
        if self.line_is_reversed:
            coll_pos = L - coll_pos

        coll_types = [line[name].__class__.__name__
                      for name in collimator_names]
        coll_weights = weights[coll_mask]
        nabs = [coll_weights[coll_losses == name].sum()
                for name in collimator_names]
        energy_weights = coll_weights * part.energy[coll_mask]
        coll_energy = [energy_weights[coll_losses == name].sum()
                       for name in collimator_names]

        return {
            's':      coll_pos,
            'name':   collimator_names,
            'n':      nabs,
            'e':      coll_energy,
            'length': coll_lengths,
            'type':   coll_types
        }

    def _add_summary(self, colldata):
        # TODO: quick lazy hack array -> list -> array. Please fix.
        data = self._summary.copy()
        data['s'] = list(data['s'])
        data['name'] = list(data['name'])
        data['n'] = list(data['n'])
        data['e'] = list(data['e'])
        data['length'] = list(data['length'])
        data['type'] = list(data['type'])
        if 'e' not in colldata:
            if not self._energy_waring_given:
                print("Warning: The JSON file does not contain energy data. "
                    + "Energy data in the loss map will be wrong.")
                self._energy_waring_given = True
            colldata['e'] = [0] * len(colldata['s'])
        if 'type' not in colldata:
            colldata['type'] = ['Unknown'] * len(colldata['s'])

        for ss, nn, nabs, ee, ll, tt in zip(colldata['s'], colldata['name'],
            colldata['n'], colldata['e'], colldata['length'], colldata['type']):
            if nn not in data['name']:
                data['s'].append(ss)
                data['name'].append(nn)
                data['n'].append(nabs)
                data['e'].append(ee)
                data['length'].append(ll)
                data['type'].append(tt)
            else:
                idx = data['name'].index(nn)
                data['n'][idx] += nabs
                data['e'][idx] += ee
                if data['type'][idx] == 'Unknown':
                    data['type'][idx] = tt
                if not np.isclose(data['length'][idx], ll):
                    raise ValueError(f"Length of {nn} is different from the one used "
                                    + "to create the loss map.")
                if data['type'][idx] != tt and tt != 'Unknown':
                    raise ValueError(f"Type of {nn} is different from the one used to "
                                    + "create the loss map.")
                if not np.isclose(data['s'][idx], ss):
                    raise ValueError(f"s position of {nn} is different from the one "
                                    + "used to create the loss map.")
        data['s'] = np.array(data['s'])
        data['name'] = np.array(data['name'])
        data['n'] = np.array(data['n'])
        data['e'] = np.array(data['e'])
        data['length'] = np.array(data['length'])
        data['type'] = np.array(data['type'])
        self._summary = data
        self._summary = self.summary.sort_values("s").to_dict('list')


    def _get_aperture_losses(self, part, line, line_shift_s, weights):
        aper_mask = part.state == 0
        if len(part.s[aper_mask]) == 0:
            return {'s': [], 'name': [], 'n': [], 'e': []}

        # Get s position per particle (lost on aperture)
        L = self.machine_length
        aper_s = np.mod(part.s[aper_mask] + line_shift_s, L)
        if self._line_is_reversed:
            aper_s = L - aper_s

        # Store names of aperture markers. Note that different s could have the same
        # name due to the interpolation.
        aper_names   = [line.element_names[i] for i in part.at_element[aper_mask]]
        name_dict    = dict(zip(aper_s, aper_names)) # TODO: not floating-point-safe and slow

        # # TODO: correctly bin (now if apertures are closer than the interpolation step,
        # # this is a bit wrong)
        # end = L - np.mod(L, self.interpolation) + self.interpolation
        # aperbins = np.linspace(0, end, np.ceil(L/self.interpolation) + 1)
        # binned = aperbins[ np.digitize(aper_s, bins=aperbins, right=False) - 1 ] + self.interpolation/2

        # Create output arrays
        aper_pos       = np.unique(aper_s)
        aper_names     = [name_dict[ss] for ss in aper_pos]
        aper_weights   = weights[aper_mask]
        aper_nabs      = [aper_weights[aper_s == s].sum() for s in aper_pos]
        energy_weights = aper_weights * part.energy[aper_mask]
        aper_energy    = [energy_weights[aper_s == s].sum() for s in aper_pos]

        return {'s': aper_pos, 'name': aper_names, 'n': aper_nabs,
                'e': aper_energy}

    def _add_aperture(self, aperdata, verbose=True):
        # TODO: quick lazy hack array -> list -> array. Please fix.
        data = self._aperture.copy()
        data['s'] = list(data['s'])
        data['name'] = list(data['name'])
        data['n'] = list(data['n'])
        data['e'] = list(data['e'])
        if 'e' not in aperdata:
            if not self._energy_waring_given:
                print("Warning: The JSON file does not contain energy data. "
                    + "Energy data in the loss map will be wrong.")
                self._energy_waring_given = True
            aperdata['e'] = [0] * len(aperdata['s'])

        for ss, nn, nabs, ee in zip(aperdata['s'], aperdata['name'],
            aperdata['n'], aperdata['e']):
            if ss not in data['s']:
                data['s'].append(ss)
                data['name'].append(nn)
                data['n'].append(nabs)
                data['e'].append(ee)
            else:
                idx = data['s'].index(ss)
                data['n'][idx] += nabs
                data['e'][idx] += ee
                if data['name'][idx] != nn:
                    raise ValueError(f"Different name {nn} at position {ss} than "
                                    + "the one used to create the loss map.")
        data['s'] = np.array(data['s'])
        data['name'] = np.array(data['name'])
        data['n'] = np.array(data['n'])
        data['e'] = np.array(data['e'])
        self._aperture = data


    @staticmethod
    def _assert_valid_json(lossmap):
        enforce_new_format = 'xcoll' in lossmap
        if 'machine_length' not in lossmap:
            raise ValueError("The JSON file does not contain the machine length data.")
        if 'interpolation' not in lossmap:
            raise ValueError("The JSON file does not contain the interpolation data.")
        if 'reversed' not in lossmap:
            raise ValueError("The JSON file does not contain the reversed data.")
        if 'momentum' not in lossmap and enforce_new_format:
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
        if 'e' not in lossmap['collimator'] and enforce_new_format:
            raise ValueError("The JSON file does not contain the collimator energy data.")
        if 'type' not in lossmap['collimator'] and enforce_new_format:
            raise ValueError("The JSON file does not contain the collimator type data.")
        if 'aperture' not in lossmap:
            raise ValueError("The JSON file does not contain the aperture data.")
        if 's' not in lossmap['aperture']:
            raise ValueError("The JSON file does not contain the aperture s data.")
        if 'name' not in lossmap['aperture']:
            raise ValueError("The JSON file does not contain the aperture name data.")
        if 'n' not in lossmap['aperture']:
            raise ValueError("The JSON file does not contain the aperture n data.")
        if 'e' not in lossmap['aperture'] and enforce_new_format:
            raise ValueError("The JSON file does not contain the aperture energy data.")


class MultiLossMap(LossMap):
    def __init__(self, *lms):
        """
        Create a MultiLossMap object from a list of LossMap objects.
        """
        self._lms = []
        self._aperture = {'s': [], 'name': [], 'n': [], 'e': []}
        self._summary = {'s': [], 'name': [], 'n': [], 'e': [], 'length': [],
                         'type': []}
        self._machine_length = None
        self._momentum = None
        self._interpolation = None
        self._xcoll = __version__
        self._energy_waring_given = False
        for lm in lms:
            self.add_lossmap(lm)


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
        self._add_summary(lm._summary)
        self._add_aperture(lm._aperture)


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

