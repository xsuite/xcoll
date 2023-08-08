# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt

from .interaction_types import source, interactions, shortcuts
from ..general import _pkg_root

import numpy as np
import pandas as pd


class CollimatorImpacts(xt.BeamElement):

    _xofields = {
        '_index':            xt.RecordIndex,
        'at_element':        xo.Int64[:],
        'at_turn':           xo.Int64[:],
        'ds':                xo.Float64[:],
        '_inter':            xo.Int64[:],
        'parent_id':         xo.Int64[:],
        'parent_x':          xo.Float64[:],
        'parent_px':         xo.Float64[:],
        'parent_y':          xo.Float64[:],
        'parent_py':         xo.Float64[:],
        'parent_zeta':       xo.Float64[:],
        'parent_delta':      xo.Float64[:],
        'parent_energy':     xo.Float64[:],
        'parent_mass':       xo.Float64[:],
        'parent_charge':     xo.Int64[:],
        'parent_z':          xo.Int64[:],
        'parent_a':          xo.Int64[:],
        'parent_pdgid':      xo.Int64[:],
        'child_id':          xo.Int64[:],
        'child_x':           xo.Float64[:],
        'child_px':          xo.Float64[:],
        'child_y':           xo.Float64[:],
        'child_py':          xo.Float64[:],
        'child_zeta':        xo.Float64[:],
        'child_delta':       xo.Float64[:],
        'child_energy':      xo.Float64[:],
        'child_mass':        xo.Float64[:],
        'child_charge':      xo.Int64[:],
        'child_z':           xo.Int64[:],
        'child_a':           xo.Int64[:],
        'child_pdgid':       xo.Int64[:],
    }

    _extra_c_sources = [
        source,
        _pkg_root.joinpath('headers','particle_states.h'),
        _pkg_root.joinpath('impacts','impacts_src','impacts.h')
    ]

    @property
    def interaction_type(self):
        return np.array([interactions[inter] for inter in self._inter])

    def collimator_name(self, element_id):
        if not hasattr(self, '_coll_ids'):
            return element_id
        elif element_id not in self._coll_ids:
            raise ValueError(f"Element {element_id} not found in list of collimators of this impact table!\n"
                            + "Did the line change without updating the list in the impact table?")
        else:
            return self._coll_ids[element_id]

    def to_pandas(self):
        n_rows = self._index.num_recorded
        coll_header = 'collimator' if hasattr(self, '_coll_ids') else 'collimator id'
        df = pd.DataFrame({
                'turn':              self.at_turn[:n_rows],
                coll_header:         [self.collimator_name(element_id) for element_id in self.at_element[:n_rows]],
                'interaction_type':  [interactions[inter] for inter in self._inter[:n_rows]],
                'ds':                self.ds[:n_rows],
                **{
                    f'{p} {val}': getattr(self, f'{p}_{val}')[:n_rows]
                    for p in ['parent', 'child']
                    for val in ['id', 'x', 'px', 'y', 'py', 'zeta', 'delta', 'energy', 'mass', 'charge', 'z', 'a', 'pdgid']
                }
            })
        return df

    # TODO: list of impacted collimators

    
    # TODO: does not work when multiple children
    # TODO: allow to select collimator by name
    def interactions_per_collimator(self, collimator_id=0, *, turn=None):
        mask = (self._inter > 0) & (self.at_element == collimator_id)
        if turn is not None:
            mask = mask & (self.at_turn == turn)
        df = pd.DataFrame({
                'int':  [shortcuts[inter] for inter in self._inter[mask]],
                'pid':               self.parent_id[mask]
            })
        return df.groupby('pid', sort=False)['int'].agg(list)
        
