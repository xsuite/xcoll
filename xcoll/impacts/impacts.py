# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt

from .interaction_types import source, interactions
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
        n_rows = self._index.num_recorded
        return np.array([interactions[inter] for inter in self._inter[:n_rows]])

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
                coll_header:         [self.collimator_name(element_id) for element_id in self.at_element[:n_rows]],
                'ds':                self.ds[:n_rows],
                'turn':              self.at_turn[:n_rows],
                'interaction_type':  self.interaction_type,
            })
        cols = ['id', 'x', 'px', 'y', 'py', 'zeta', 'delta', 'energy', 'mass', 'charge', 'z', 'a', 'pdgid']
        for particle in ['parent', 'child']:
            multicols = pd.MultiIndex.from_tuples([(particle, col) for col in cols])
            newdf = pd.DataFrame(index=df.index, columns=multicols)
            for col in cols:
                newdf[particle, col] = getattr(self, particle + '_' + col)[:n_rows]
            df = pd.concat([df, newdf], axis=1)
        return df
