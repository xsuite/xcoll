# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt

from .interaction_types import source, interactions, shortcuts, is_point
from ..general import _pkg_root

import numpy as np
import pandas as pd


class InteractionRecord(xt.BeamElement):
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

    allow_track = False

    _extra_c_sources = [
        source,
        _pkg_root.joinpath('headers','particle_states.h'),
        _pkg_root.joinpath('interaction_record','interaction_record_src','interaction_record.h')
    ]


    @classmethod
    def start(cls, line, names=None, *, record_touches=None, record_scatterings=None, capacity=1e6, io_buffer=None):
        names = _get_xcoll_elements(line, names)
        if len(names) == 0:
            return
        elements = [line[name] for name in names]
        capacity = int(capacity)
        if io_buffer is None:
            io_buffer = xt.new_io_buffer(capacity=capacity)
        if record_touches is None and record_scatterings is None:
            record_touches = True
            record_scatterings = True
        elif record_touches is None:
            record_touches = not record_scatterings
        elif record_scatterings is None:
            record_scatterings = not record_touches
        assert record_touches is True or record_touches is False
        assert record_scatterings is True or record_scatterings is False
        for el in elements:
            if not el.record_touches and not el.record_scatterings:
                el.record_touches = record_touches
                el.record_scatterings = record_scatterings
        record = xt.start_internal_logging(io_buffer=io_buffer, capacity=capacity, \
                                           elements=elements)
        record._line = line
        record._io_buffer = io_buffer
        record._recording_elements = names
        record._coll_ids = {name: line.element_names.index(name) for name in names}
        record._coll_names = {vv: kk for kk, vv in record._coll_ids.items()}
        return record

    def stop(self, names=None):
        self.assert_class_init()
        names = _get_xcoll_elements(self.line, names)
        elements = [self.line[name] for name in names]
        if self.line.tracker is not None:
            self.line.tracker._check_invalidated()
        xt.stop_internal_logging(elements=elements)
        # Removed the stopped collimators from list of logged elements
        self._recording_elements = list(set(self._recording_elements) - set(names))


    def assert_class_init(self):
        if not hasattr(self, '_io_buffer') or not hasattr(self, '_line') \
        or not hasattr(self, '_recording_elements'):
            raise ValueError("This InteractionRecord has been manually instantiated, "
                           + "hence the expanded API is not available. Use "
                           + "InteractionRecord.start() to initialise with extended API.")

    @property
    def line(self):
        if hasattr(self, '_line'):
            return self._line

    @property
    def io_buffer(self):
        if hasattr(self, '_io_buffer'):
            return self._io_buffer

    @property
    def capacity(self):
        if hasattr(self, '_io_buffer'):
            return self.io_buffer.capacity

    # @capacity.setter
    # def capacity(self, val):
    #     if hasattr(self, '_io_buffer'):
    #         capacity = int(capacity)
    #         if capacity < self.capacity:
    #             raise NotImplementedError("Shrinking of capacity not yet implemented!")
    #         elif capacity == self.capacity:
    #             return
    #         else:
    #             self.io_buffer.grow(capacity - self.capacity)
    #             # TODO: increase capacity of iobuffer AND of fields in record table

    @property
    def recording_elements(self):
        if hasattr(self, '_recording_elements'):
            return self._recording_elements

    @recording_elements.setter
    def recording_elements(self, val):
        self.assert_class_init()
        if val is None:
            val = []
        record_start = _get_xcoll_elements(self.line, val)
        self.stop(set(self.recording_elements) - set(record_start))
        elements = [self.line[name] for name in record_start]
        for el in elements:
            if not el.record_touches and not el.record_scatterings:
                el.record_touches = True
                el.record_scatterings = True
        xt.start_internal_logging(io_buffer=self.io_buffer, capacity=self.capacity, \
                                  record=self, elements=elements)
        self._recording_elements = record_start
        # Updating coll IDs: careful to correctly overwrite existing values
        self._coll_ids.update({name: self.line.element_names.index(name) for name in record_start})
        self._coll_names = {vv: kk for kk, vv in self._coll_ids.items()}

    @property
    def interaction_type(self):
        return np.array([interactions[inter] for inter in self._inter])

    def _collimator_name(self, element_id):
        if not hasattr(self, '_coll_names'):
            return element_id
        elif element_id not in self._coll_names:
            raise ValueError(f"Element {element_id} not found in list of collimators of this record table! "
                           + f"Did the line change without updating the list in the table?")
        else:
            return self._coll_names[element_id]

    def _collimator_id(self, element_name):
        if not hasattr(self, '_coll_ids'):
            return element_name
        elif element_name not in self._coll_ids:
            raise ValueError(f"Element {element_name} not found in list of collimators of this record table! "
                           + f"Did the line change without updating the list in the table?")
        else:
            return self._coll_ids[element_name]

    def to_pandas(self):
        n_rows = self._index.num_recorded
        coll_header = 'collimator' if hasattr(self, '_coll_names') else 'collimator_id'
        df = pd.DataFrame({
                'turn':              self.at_turn[:n_rows],
                coll_header:         [self._collimator_name(element_id) for element_id in self.at_element[:n_rows]],
                'interaction_type':  [interactions[inter] for inter in self._inter[:n_rows]],
                'ds':                self.ds[:n_rows],
                **{
                    f'{p}_{val}': getattr(self, f'{p}_{val}')[:n_rows]
                    for p in ['parent', 'child']
                    for val in ['id', 'x', 'px', 'y', 'py', 'zeta', 'delta', 'energy', 'mass', 'charge', 'z', 'a', 'pdgid']
                }
            })
        return df

    # TODO: list of impacted collimators


    # TODO: does not work when multiple children
    def interactions_per_collimator(self, collimator=0, *, turn=None):
        if isinstance(collimator, str):
            collimator = self._collimator_id(collimator)
        mask = (self._inter > 0) & (self.at_element == collimator)
        if turn is not None:
            mask = mask & (self.at_turn == turn)
            df = pd.DataFrame({
                    'int':  [shortcuts[inter] for inter in self._inter[mask]],
                    'pid':  self.parent_id[mask]
                })
            return df.groupby('pid', sort=False)['int'].agg(list)
        else:
            df = pd.DataFrame({
                    'int':   [shortcuts[inter] for inter in self._inter[mask]],
                    'turn':  self.at_turn[mask],
                    'pid':   self.parent_id[mask]
                })
            return df.groupby(['pid', 'turn'], sort=False)['int'].apply(list)

    def first_touch_per_turn(self):
        n_rows = self._index.num_recorded
        df = pd.DataFrame({'parent_id': self.parent_id[:n_rows],
                           'at_turn': self.at_turn[:n_rows],
                           'at_element': self.at_element[:n_rows]})
        mask = np.char.startswith(self.interaction_type[:n_rows], 'Enter Jaw')
        idx_first = [group.at_element.idxmin() for _, group in df[mask].groupby(['at_turn', 'parent_id'], sort=False)]
        df_first = self.to_pandas().loc[idx_first]
        df_first.insert(2, "jaw", df_first.interaction_type.astype(str).str[-1])
        to_drop = ['ds', 'interaction_type',
                   *[col for col in df_first.columns if col.startswith('child_')]]
        to_rename = {col: col.replace('parent_', '') for col in df_first.columns if col.startswith('parent_')}
        return df_first.drop(columns=to_drop).rename(columns=to_rename)


def _get_xcoll_elements(line, names):
    from xcoll import element_classes
    if names is None or names is True:
        names = line.get_elements_of_type(element_classes)[1]
        if len(names) == 0:
            raise ValueError("No Xcoll elements in line!")
    if names is False:
        names = []
    if not hasattr(names, '__iter__') or isinstance(names, str):
        names = [names]
    for name in names:
        if name not in line.element_names:
            raise ValueError(f"Element {name} not found in line!")
        if not isinstance(line[name], element_classes):
            raise ValueError(f"Element {name} not an Xcoll element!")
    return names

