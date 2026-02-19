# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import xobjects as xo
import xtrack as xt

from .interaction_types import interactions_src, interaction_names, shortcuts
from ..general import _pkg_root
from ..headers.particle_states import particle_states_src

import numpy as np
import pandas as pd

interaction_names = {kk: vv.replace('_', ' ').title().\
                            replace('Pn ','PN ').replace('Pp ','PP ').replace(' Mcs',' MCS').\
                            replace(' Ch',' CH').replace(' Vr',' VR')
                     for kk, vv in interaction_names.items()}


class InteractionRecord(xt.BeamElement):
    _xofields = {
        '_index':            xt.RecordIndex,
        'at_element':        xo.Int64[:],
        'at_turn':           xo.Int64[:],
        '_inter':            xo.Int64[:],
        'id_before':         xo.Int64[:],
        's_before':          xo.Float64[:],
        'x_before':          xo.Float64[:],
        'px_before':         xo.Float64[:],
        'y_before':          xo.Float64[:],
        'py_before':         xo.Float64[:],
        'zeta_before':       xo.Float64[:],
        'delta_before':      xo.Float64[:],
        'energy_before':     xo.Float64[:],
        'mass_before':       xo.Float64[:],
        'charge_before':     xo.Int64[:],
        'z_before':          xo.Int64[:],
        'a_before':          xo.Int64[:],
        'pdgid_before':      xo.Int64[:],
        'id_after':          xo.Int64[:],
        's_after':           xo.Float64[:],
        'x_after':           xo.Float64[:],
        'px_after':          xo.Float64[:],
        'y_after':           xo.Float64[:],
        'py_after':          xo.Float64[:],
        'zeta_after':        xo.Float64[:],
        'delta_after':       xo.Float64[:],
        'energy_after':      xo.Float64[:],
        'mass_after':        xo.Float64[:],
        'charge_after':      xo.Int64[:],
        'z_after':           xo.Int64[:],
        'a_after':           xo.Int64[:],
        'pdgid_after':       xo.Int64[:],
    }

    allow_track = False

    _extra_c_sources = [
        interactions_src,
        particle_states_src,
        _pkg_root.joinpath('interaction_record','interaction_record_src','interaction_record.h')
    ]


    @classmethod
    def start(cls, *, line=None, elements=None, names=None, record_impacts=None, record_exits=None,
              record_scatterings=None, capacity=1e6, io_buffer=None, coll_ids=None):
        elements, names = _get_xcoll_elements(line, elements, names)
        if len(names) == 0:
            return
        capacity = int(capacity)

        if getattr(line, 'tracker', None) is None \
        or getattr(line.tracker, 'io_buffer', None) is None:
            if io_buffer is None:
                io_buffer = xt.new_io_buffer(capacity=capacity)
        elif io_buffer is not None:
            raise ValueError("Cannot provide io_buffer when tracker already built!")
        else:
            io_buffer = line.tracker.io_buffer
        if capacity > io_buffer.capacity:
            io_buffer.grow(capacity - io_buffer.capacity)
        if record_impacts is None and record_scatterings is None:
            record_impacts = True
            record_scatterings = True
        elif record_impacts is None:
            record_impacts = not record_scatterings
        elif record_scatterings is None:
            record_scatterings = not record_impacts
        if record_exits is None:
            # record_exits defaults to True only if the other two are True
            record_exits = record_impacts and record_scatterings
        assert record_impacts is True or record_impacts is False
        assert record_exits is True or record_exits is False
        assert record_scatterings is True or record_scatterings is False
        for el in elements:
            if not el.record_impacts and not el.record_exits and not el.record_scatterings:
                el.record_impacts = record_impacts
                el.record_exits = record_exits
                el.record_scatterings = record_scatterings
        record = xt.start_internal_logging(io_buffer=io_buffer, capacity=capacity, \
                                           elements=elements)
        record._line = line
        record._io_buffer = io_buffer
        recording_elements = names if len(names) > 0 else elements
        record._recording_elements = recording_elements
        if coll_ids is None:
            if line is None:
                if len(names) > 0:
                    record._coll_ids = {name: idx for idx, name in enumerate(names)}
            else:
                record._coll_ids = {name: line.element_names.index(name) for name in names}
        else:
            assert len(coll_ids) == len(names)
            record._coll_ids = {name: idx for name, idx in zip(names, coll_ids)}
        if hasattr(record, '_coll_ids'):
            record._coll_names = {vv: kk for kk, vv in record._coll_ids.items()}
        return record

    def stop(self, *, elements=False, names=False):
        self.assert_class_init()
        elements, names = _get_xcoll_elements(self.line, elements, names)
        if self.line is not None and self.line.tracker is not None:
            self.line.tracker._check_invalidated()
        xt.stop_internal_logging(elements=elements)
        # Removed the stopped collimators from list of logged elements
        stopping_elements = names if len(names) > 0 else elements
        self._recording_elements = list(set(self._recording_elements) - set(stopping_elements))


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
        for el in elements:  # TODO: this should be smarter
            if not el.record_impacts and not el.record_scatterings:
                el.record_impacts = True
                el.record_scatterings = True
        xt.start_internal_logging(io_buffer=self.io_buffer, capacity=self.capacity, \
                                  record=self, elements=elements)
        self._recording_elements = record_start
        # Updating coll IDs: careful to correctly overwrite existing values
        self._coll_ids.update({name: self.line.element_names.index(name) for name in record_start})
        self._coll_names = {vv: kk for kk, vv in self._coll_ids.items()}

    @property
    def interaction_type(self):
        return np.array([interaction_names[inter] for inter in self._inter])

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

    def to_pandas(self, frame=None):
        if frame is None:
            frame = 'jaw'
        frame = frame.lower()
        if frame not in ['jaw', 'collimator', 'lab']:
            raise ValueError(f"Invalid frame {frame}. Must be 'jaw', 'collimator', or 'lab'!")
        n_rows = self._index.num_recorded
        coll_header = 'collimator' if hasattr(self, '_coll_names') else 'collimator_id'
        df = pd.DataFrame({
                'turn':              self.at_turn[:n_rows],
                coll_header:         [self._collimator_name(element_id) for element_id in self.at_element[:n_rows]],
                'interaction_type':  [interaction_names[inter] for inter in self._inter[:n_rows]]
            } | {
                f'{val}_{p}': getattr(self, f'{val}_{p}')[:n_rows]
                for p in ['before', 'after']
                for val in ['id', 's', 'x', 'px', 'y', 'py', 'zeta', 'delta', 'energy', 'mass', 'charge', 'z', 'a', 'pdgid']
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
                    'pid':  self.id_before[mask]
                })
            return df.groupby('pid', sort=False, group_keys=False)['int'].agg(list)
        else:
            df = pd.DataFrame({
                    'int':   [shortcuts[inter] for inter in self._inter[mask]],
                    'turn':  self.at_turn[mask],
                    'pid':   self.id_before[mask]
                })
            return df.groupby(['pid', 'turn'], sort=False, group_keys=False)['int'].apply(list)

    def first_touch_per_turn(self, frame=None):
        n_rows = self._index.num_recorded
        df = pd.DataFrame({'id_before': self.id_before[:n_rows],
                           'at_turn': self.at_turn[:n_rows],
                           'at_element': self.at_element[:n_rows]})
        mask = np.char.startswith(self.interaction_type[:n_rows], 'Enter Jaw')
        idx_first = [group.at_element.idxmin() for _, group in df[mask].groupby(
                        ['at_turn', 'id_before'], sort=False, group_keys=False)]
        df_first = self.to_pandas(frame=frame).loc[idx_first]
        df_first.insert(2, "jaw", df_first.interaction_type.astype(str).str[-1])
        to_drop = ['interaction_type',
                   *[col for col in df_first.columns if col.endswith('_after')]]
        to_rename = {col: col.replace('_before', '') for col in df_first.columns if col.endswith('before')}
        to_rename['id_before'] = 'pid'
        return df_first.drop(columns=to_drop).rename(columns=to_rename)


def _get_xcoll_elements(line=None, elements=None, names=None):
    from xcoll.beam_elements import block_classes
    if names is not None and names is not False and \
    (not hasattr(names, '__iter__') or isinstance(names, str)):
        names = [names]
    if elements is not None and elements is not False and \
    (not hasattr(elements, '__iter__') or isinstance(elements, str)):
        elements = [elements]
    if line is None:
        if elements is None:
            raise ValueError("No line nor elements provided!")
    else:
        if elements is not None and elements is not False:
            raise ValueError("Cannot provide both line and elements!")
        if names is None or names is True:
            elements, names = line.get_elements_of_type(block_classes)
            if len(names) == 0:
                raise ValueError("No Xcoll elements in line!")
        elif names is False:
            names = []
            elements = []
        else:
            assert elements is not False
            for name in names:
                if name not in line.element_names:
                    raise ValueError(f"Element {name} not found in line!")
            elements = [line[name] for name in names]
    for idx, element in enumerate(elements):
        if not isinstance(element, block_classes):
            name = name[idx] if names is not None else element.__class__.__name__
            raise ValueError(f"Element {name} not an Xcoll element!")
    return elements, names
