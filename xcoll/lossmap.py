# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pandas as pd
from pathlib import Path
import json

import xtrack as xt
import xpart as xp
import xobjects as xo

from .beam_elements import collimator_classes


class LossMap:

    def __init__(self, line, part, *, line_is_reversed, interpolation=0.1,
                 weights=None, weight_function=None):
        self._line = line
        self._line_is_reversed = line_is_reversed
        self._machine_length = line.get_length()
        self._part = part
        self._interpolation = interpolation
        if weights is None:
            if weight_function is None:
                self._weights = np.ones(len(part.x))
            else:
                self._weights = _create_weights_from_initial_state(part, weight_function)
        else:
            if weight_function is not None:
                raise ValueError("Use either 'weights' or 'weight_function', not both!")
            self._weights = part.sort(interleave_lost_particles=True)

        # Correct particles that are lost in aperture directly after collimator -> should be absorbed
        self._correct_absorbed()

        # Loss location refinement
        if interpolation is not None:
            self._interpolate()

        self._make_coll_summary()
        coll_summary = self._summary[self._summary.nabs > 0].to_dict('list')
        aper_s, aper_names, aper_nabs, aper_energy = self._get_aperture_losses()

        self._lossmap = {
                'collimator': {
                    's':      coll_summary['s'],
                    'name':   coll_summary['collname'],
                    'length': coll_summary['length'],
                    'n':      coll_summary['nabs'],
                    # 'e':      coll_summary['energy']
                }
                ,
                'aperture': {
                    's':    aper_s,
                    'name': aper_names,
                    'n':    aper_nabs,
                    # 'e':    aper_energy
                }
                ,
                'machine_length': self._machine_length
                ,
                'interpolation': interpolation
                ,
                'reversed': self._line_is_reversed
            }


    def to_json(self, file):
        with open(Path(file), 'w') as fid:
            json.dump(self._lossmap, fid, indent=True, cls=xo.JEncoder)

    def save_summary(self, file):
        with open(Path(file), 'w') as fid:
            fid.write(self._summary.__repr__())


    @property
    def lossmap(self):
        return self._lossmap

    @property
    def summary(self):
        return self._summary

    @property
    def line(self):
        return self._line

    @property
    def line_is_reversed(self):
        return self._line_is_reversed

    @property
    def machine_length(self):
        return self._machine_length

    @property
    def part(self):
        return self._part

    @property
    def interpolation(self):
        return self._interpolation

    @property
    def weights(self):
        return self._weights


    def _correct_absorbed(self):
        # Correct particles that are at an aperture directly after a collimator
        part = self._part
        for idx, elem in enumerate(part.at_element):
            if (part.state[idx] == 0 and elem > 0
            and self._line.element_names[elem-1] in
            self._line.get_elements_of_type(collimator_classes)[1]):
                print(f"Found at {self._line.element_names[elem]}, "
                    + f"should be {self._line.element_names[elem-1]}")
                part.at_element[idx] = elem - 1
                what_type = self._line[elem-1].__class__.__name__
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

    def _interpolate(self):
        aper_s = list(self._part.s[self._part.state==0])
        if len(aper_s) > 0:
            print("Performing the aperture losses refinement.")
            loss_loc_refinement = xt.LossLocationRefinement(
                self._line,
                n_theta = 360,            # Angular resolution
                r_max = 0.5,              # Maximum transverse aperture [m]
                dr = 50e-6,               # Transverse accuracy [m]
                ds = self.interpolation   # Longitudinal accuracy [m]
            )
            loss_loc_refinement.refine_loss_location(self._part)


    def _make_coll_summary(self):
        collimator_names = self._line.get_elements_of_type(collimator_classes)[1]
        coll_mask     = (self._part.state <= -330) & (self._part.state >= -340)
        coll_losses   = np.array([self._line.element_names[i]
                                  for i in self._part.at_element[coll_mask]])
        coll_lengths  = [self._line[j].length for j in collimator_names] 
        coll_pos      = [(self._line.get_s_position(i) + self._line[i].length/2)
                         for i in collimator_names]

        if self._line is reversed:
            coll_pos  = [self._machine_length - s for s in coll_pos]

        coll_types    = [self._line[i].__class__.__name__ for i in collimator_names]  
        coll_weights  = self._weights[coll_mask]
        nabs          = [coll_weights[coll_losses == j].sum() for j in collimator_names]

        self._summary = pd.DataFrame({
            'collname': collimator_names,
            'nabs':     nabs,    # of particles lost on collimators
            # 'energy':   energy,
            'length':   coll_lengths,
            's':        coll_pos,
            'type':     coll_types
        })


    def _get_aperture_losses(self):
        # Get s position per particle (lost on aperture)
        aper_mask    = self._part.state == 0
        aper_s       = list(self._part.s[aper_mask])

        if len(aper_s) == 0:
            return [], [], []
        if self._line_is_reversed:
            aper_s   = [ self._machine_length - s for s in aper_s ]

        # Store names of aperture markers
        aper_names   = [self._line.element_names[i] for i in self._part.at_element[aper_mask]]
        name_dict    = dict(zip(aper_s, aper_names)) # TODO: not floating-point-safe and slow

        # Create output arrays
        aper_pos     = np.unique(aper_s)
        aper_weights = self._weights[aper_mask]
        aper_nabs    = [aper_weights[aper_s == j].sum() for j in aper_pos] 
        aper_names   = [name_dict[ss] for ss in aper_pos]

        aper_energy = 0

        return aper_pos, aper_names, aper_nabs, aper_energy


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

