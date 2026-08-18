# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import pandas as pd
from warnings import warn

import xobjects as xo
import xtrack as xt

from .interaction_types import interactions_src, interaction_names, shortcuts
from ..general import _pkg_root
from ..headers.particle_states import particle_states_src

interaction_names = {kk: vv.replace('_', ' ').title().\
                            replace('Pn ','PN ').replace('Pp ','PP ').\
                            replace(' Mcs',' MCS').replace(' Ch',' CH').\
                            replace(' Vr',' VR')
                     for kk, vv in interaction_names.items()}


class InteractionRecord(xt.BeamElement):
    _xofields = {
        '_index':              xt.RecordIndex,
        '_record_all_columns': xo.Int8,
        'at_turn':             xo.Int64[:],
        'at_element':          xo.Int64[:],
        'shape_id':            xo.Int64[:],
        '_inter':              xo.Int64[:],
        'particle_id_before':  xo.Int64[:],
        's_before':            xo.Float64[:],
        'x_before':            xo.Float64[:],
        'px_before':           xo.Float64[:],
        'y_before':            xo.Float64[:],
        'py_before':           xo.Float64[:],
        'zeta_before':         xo.Float64[:],
        'delta_before':        xo.Float64[:],
        'energy_before':       xo.Float64[:],
        'mass_before':         xo.Float64[:],
        'pdgid_before':        xo.Int64[:],
        'particle_id_after':   xo.Int64[:],
        's_after':             xo.Float64[:],
        'x_after':             xo.Float64[:],
        'px_after':            xo.Float64[:],
        'y_after':             xo.Float64[:],
        'py_after':            xo.Float64[:],
        'zeta_after':          xo.Float64[:],
        'delta_after':         xo.Float64[:],
        'energy_after':        xo.Float64[:],
        'mass_after':          xo.Float64[:],
        'pdgid_after':         xo.Int64[:],
    }

    allow_track = False

    _extra_c_sources = [
        interactions_src,
        particle_states_src,
        _pkg_root.joinpath('interaction_record','interaction_record_src','interaction_record.h')
    ]

    def __init__(self, *, line=None, elements=None, names=None, num_rows=1e6,
                 columns=None, record_impacts=None, record_exits=None,
                 record_scatterings=None, io_buffer=None, **kwargs):
        elements, names = _get_xcoll_elements(line, elements, names)
        if len(elements) == 0:
            raise ValueError("No Xcoll elements provided to InteractionRecord!")

        # Initialise record table
        if 'capacity' in kwargs:
            warn("`capacity` in the InteractionRecord init is deprecated. Use "
                 "`num_rows` instead.", FutureWarning)
            num_rows = kwargs.pop('capacity')
        num_rows = int(num_rows)
        init_dict, record_all_columns, recorded_columns = _make_columns(
            num_rows=num_rows, columns=columns
        )
        init_dict.update(kwargs)

        # Memory layout of InteractionRecord:
        #       0: size
        #       0: _index (size 32)
        #      32: _record_all_columns (size 8)
        #   40-96: <8 transformations (shift_x, rot_x_rad, ...)> (size 8*8=64)
        # 104-304: <pointers to all 26 arrays> (size 26*8=208)
        # 312-...: <all arrays> (each has size 8*num_rows + 16)
        #
        # Total size: 312 + num_logged*8*num_rows + 26*16

        # Get context and buffer with correct capacity
        _context = kwargs.get('_context')
        capacity = 728 + 8*num_rows*len(recorded_columns)
        if getattr(line, 'tracker', None) is None \
        or getattr(line.tracker, 'io_buffer', None) is None:
            if io_buffer is None:
                io_buffer = xt.new_io_buffer(_context=_context,
                                             capacity=capacity)
        elif io_buffer is not None:
            raise ValueError("Cannot provide io_buffer when tracker already built!")
        else:
            io_buffer = line.tracker.io_buffer
            if _context is None:
                init_dict['_context'] = line.tracker._context
        if io_buffer.capacity < capacity:
            io_buffer.grow(capacity - io_buffer.capacity)
        if _context is not None and type(_context) is not \
        type(io_buffer.context):
                raise ValueError("io_buffer context does not match provided "
                                 "context!")

        super().__init__(_buffer=io_buffer, **init_dict)

        _set_recording_flags(elements, record_impacts, record_exits,
                             record_scatterings)

        self._line = line
        self._io_buffer = io_buffer
        self._index.capacity = num_rows
        self._record_all_columns = record_all_columns
        self._recorded_columns = recorded_columns
        self._recording_elements = {
            name: el for name, el in zip(names, elements)
        }
        if line is None:
            self._coll_ids = {name: idx for idx, name in enumerate(names)}
        else:
            self._coll_ids = {
                name: line.element_names.index(name) for name in names
            }
        self._coll_names = {vv: kk for kk, vv in self._coll_ids.items()}

        xt.start_internal_logging(io_buffer=io_buffer, record=self,
                                  elements=elements)


    @classmethod
    def start(cls, *, line=None, elements=None, names=None, record_impacts=None,
              record_exits=None, record_scatterings=None, capacity=1e6,
              io_buffer=None, columns=None):
        warn("InteractionRecord.start() is deprecated. Use "
             "InteractionRecord(...) instead.", FutureWarning)
        return cls(line=line, elements=elements, names=names, capacity=capacity,
                   record_impacts=record_impacts, record_exits=record_exits,
                   record_scatterings=record_scatterings, io_buffer=io_buffer,
                   columns=columns)

    def stop(self, *, elements=None, names=None):
        elements, names = _get_xcoll_elements(self.line, elements, names)
        if self.line is not None and self.line.tracker is not None:
            self.line.tracker._check_invalidated()
        xt.stop_internal_logging(elements=elements)
        # Removed the stopped collimators from list of logged elements
        stopping_elements = elements
        self._recording_elements = {
            name: el for name, el in self._recording_elements.items()
            if el not in stopping_elements
        }

    @property
    def size(self):
        return self._xobject._size

    @property
    def line(self):
        if hasattr(self, '_line'):
            return self._line

    @property
    def io_buffer(self):
        if hasattr(self, '_io_buffer'):
            return self._io_buffer

    @property
    def num_rows(self):
        return self._index.capacity

    @property
    def io_buffer_capacity(self):
        if hasattr(self, '_io_buffer'):
            return self.io_buffer.capacity

    @property
    def recording_elements(self):
        if hasattr(self, '_recording_elements'):
            return self._recording_elements

    @property
    def interaction_type(self):
        return np.array([interaction_names[inter] for inter in self._inter])

    def to_pandas(self, frame=None):
        if frame is None:
            frame = 'jaw'
        frame = frame.lower()
        if frame not in ['jaw', 'collimator', 'lattice']:
            raise ValueError(f"Invalid frame {frame}. Must be 'jaw', "
                             f"'collimator', or 'lattice'!")
        n_rows = self._index.num_recorded
        nparr = self._context.nparray_from_context_array
        interaction_type = [
                    inter.tolist() if hasattr(inter, 'tolist') else inter
                    for inter in self._inter[:n_rows]
        ]
        data = {
            'turn':             nparr(self.at_turn[:n_rows]),
            'collimator':       [self._collimator_name(element_id)
                                 for element_id in self.at_element[:n_rows]],
            'interaction_type': [interaction_names[inter]
                                 for inter in interaction_type],
        }
        for p in ['before', 'after']:
            for val in ['particle_id', 's', 'x', 'px', 'y', 'py', 'zeta',
                        'delta', 'energy', 'mass', 'pdgid']:
                field = f'{val}_{p}'
                if self._column_is_recorded(field):
                    data[field] = nparr(getattr(self, field)[:n_rows])
        df = pd.DataFrame(data)

        # Different reference frames:
        #   - lattice:       as it is in the lattice (pipe frame)
        #   - collimator:    rotated to the (potentially jaw-dependent)
        #                    collimator angle
        #   - jaw (default): as collimator, but also rotated to the tilt
        #                    angle, moved to the upstream jaw corner,
        #                    and mirrored for the right jaw
        if frame != 'jaw':
            required_columns = [
                f'{val}_{p}' for p in ['before', 'after']
                for val in ['s', 'x', 'px', 'delta']
            ]
            if frame == 'lattice':
                required_columns += [
                    f'{val}_{p}' for p in ['before', 'after']
                    for val in ['y', 'py']
                ]
            self._check_columns_recorded(required_columns, frame)

            # Move back to the collimator frame

            # Coordinate arrays
            s_before     = df["s_before"].to_numpy(copy=True)
            x_before     = df["x_before"].to_numpy(copy=True)
            px_before    = df["px_before"].to_numpy(copy=True)
            y_before     = df["y_before"].to_numpy(copy=True)
            py_before    = df["py_before"].to_numpy(copy=True)
            delta_before = df["delta_before"].to_numpy(copy=False)  # No need to write to this array, so no need to copy
            s_after      = df["s_after"].to_numpy(copy=True)
            x_after      = df["x_after"].to_numpy(copy=True)
            px_after     = df["px_after"].to_numpy(copy=True)
            y_after      = df["y_after"].to_numpy(copy=True)
            py_after     = df["py_after"].to_numpy(copy=True)
            delta_after  = df["delta_after"].to_numpy(copy=False)   # No need to write to this array, so no need to copy

            # Collimator attribute arrays
            cat = df['collimator'].astype("category")
            codes = cat.cat.codes.to_numpy(copy=False)
            names = cat.cat.categories
            sh  = self.shape_id[:n_rows]
            els = self.recording_elements
            sin_zL = np.array([els[name]._sin_zL for name in names])
            cos_zL = np.array([els[name]._cos_zL for name in names])
            sin_zR = np.array([els[name]._sin_zR for name in names])
            cos_zR = np.array([els[name]._cos_zR for name in names])
            sin_yL = np.array([els[name]._sin_yL for name in names])
            cos_yL = np.array([els[name]._cos_yL for name in names])
            tilt_L = np.array([els[name].tilt_L for name in names])
            sin_yR = np.array([els[name]._sin_yR for name in names])
            cos_yR = np.array([els[name]._cos_yR for name in names])
            tilt_R = np.array([els[name].tilt_R for name in names])
            length = np.array([els[name].length for name in names])
            jaw_L  = np.array([els[name].jaw_LU for name in names])
            jaw_R  = np.array([els[name].jaw_RU for name in names])
            length_front = np.array([
                                els[name].length_front
                                if hasattr(els[name], 'length_front') else 0
                                for name in names
                            ])

            # Mirror back if on a right jaw (negative shape_id)
            idx = np.flatnonzero((x_before != -1) & (sh < 0))
            x_before[idx] = -x_before[idx]
            idx = np.flatnonzero((px_before != -1) & (sh < 0))
            px_before[idx] = -px_before[idx]    # Also valid for exact angle
            idx = np.flatnonzero((x_after != -1) & (sh < 0))
            x_after[idx] = -x_after[idx]
            idx = np.flatnonzero((px_after != -1) & (sh < 0))
            px_after[idx] = -px_after[idx]

            # Rotate back from the tilted frame
            idx = np.flatnonzero((s_before != -1) & (x_before != -1) & (sh >= 0))
            new_s = s_before[idx]*cos_yL[codes[idx]] - x_before[idx]*sin_yL[codes[idx]]
            new_x = s_before[idx]*sin_yL[codes[idx]] + x_before[idx]*cos_yL[codes[idx]]
            s_before[idx] = new_s
            x_before[idx] = new_x
            idx = np.flatnonzero((px_before != -1) & (sh >= 0))
            px_before[idx] = px_before[idx] + tilt_L[codes[idx]]*(1 + delta_before[idx])
            # xp = px/sqrt((1+delta)**2 - px*px - py*py)
            # px = xp*(1+delta)/sqrt(1 + xp*xp + yp*yp)
            idx = np.flatnonzero((s_before != -1) & (x_before != -1) & (sh < 0))
            new_s = s_before[idx]*cos_yR[codes[idx]] - x_before[idx]*sin_yR[codes[idx]]
            new_x = s_before[idx]*sin_yR[codes[idx]] + x_before[idx]*cos_yR[codes[idx]]
            s_before[idx] = new_s
            x_before[idx] = new_x
            idx = np.flatnonzero((px_before != -1) & (sh < 0))
            px_before[idx] = px_before[idx] + tilt_R[codes[idx]]*(1 + delta_before[idx])
            # xp = px/sqrt((1+delta)**2 - px*px - py*py)
            # px = xp*(1+delta)/sqrt(1 + xp*xp + yp*yp)
            idx = np.flatnonzero((s_after != -1) & (x_after != -1) & (sh >= 0))
            new_s = s_after[idx]*cos_yL[codes[idx]] - x_after[idx]*sin_yL[codes[idx]]
            new_x = s_after[idx]*sin_yL[codes[idx]] + x_after[idx]*cos_yL[codes[idx]]
            s_after[idx] = new_s
            x_after[idx] = new_x
            idx = np.flatnonzero((px_after != -1) & (sh >= 0))
            px_after[idx] = px_after[idx] + tilt_L[codes[idx]]*(1 + delta_after[idx])
            # xp = px/sqrt((1+delta)**2 - px*px - py*py)
            # px = xp*(1+delta)/sqrt(1 + xp*xp + yp*yp)
            idx = np.flatnonzero((s_after != -1) & (x_after != -1) & (sh < 0))
            new_s = s_after[idx]*cos_yR[codes[idx]] - x_after[idx]*sin_yR[codes[idx]]
            new_x = s_after[idx]*sin_yR[codes[idx]] + x_after[idx]*cos_yR[codes[idx]]
            s_after[idx] = new_s
            x_after[idx] = new_x
            idx = np.flatnonzero((px_after != -1) & (sh < 0))
            px_after[idx] = px_after[idx] + tilt_R[codes[idx]]*(1 + delta_after[idx])
            # xp = px/sqrt((1+delta)**2 - px*px - py*py)
            # px = xp*(1+delta)/sqrt(1 + xp*xp + yp*yp)

            # Move back from the jaw corner
            idx = np.flatnonzero((x_before != -1) & (sh >= 0))
            x_before[idx] += jaw_L[codes[idx]]
            idx = np.flatnonzero((s_before != -1) & (sh >= 0))
            s_before[idx] -= length[codes[idx]]/2*(1 - cos_yL[codes[idx]])
            idx = np.flatnonzero((x_before != -1) & (sh < 0))
            x_before[idx] += jaw_R[codes[idx]]
            idx = np.flatnonzero((s_before != -1) & (sh < 0))
            s_before[idx] -= length[codes[idx]]/2*(1 - cos_yR[codes[idx]])
            idx = np.flatnonzero((x_after != -1) & (sh >= 0))
            x_after[idx] += jaw_L[codes[idx]]
            idx = np.flatnonzero((s_after != -1) & (sh >= 0))
            s_after[idx] -= length[codes[idx]]/2*(1 - cos_yL[codes[idx]])
            idx = np.flatnonzero((x_after != -1) & (sh < 0))
            x_after[idx] += jaw_R[codes[idx]]
            idx = np.flatnonzero((s_after != -1) & (sh < 0))
            s_after[idx] -= length[codes[idx]]/2*(1 - cos_yR[codes[idx]])

        if frame == 'lattice':
            # Rotate back from the collimator frame
            idx = np.flatnonzero((x_before != -1) & (y_before != -1) & (sh >= 0))
            new_x = x_before[idx]*cos_zL[codes[idx]] - y_before[idx]*sin_zL[codes[idx]]
            new_y = x_before[idx]*sin_zL[codes[idx]] + y_before[idx]*cos_zL[codes[idx]]
            x_before[idx] = new_x
            y_before[idx] = new_y
            idx = np.flatnonzero((px_before != -1) & (py_before != -1) & (sh >= 0))
            new_px = px_before[idx]*cos_zL[codes[idx]] - py_before[idx]*sin_zL[codes[idx]]
            new_py = px_before[idx]*sin_zL[codes[idx]] + py_before[idx]*cos_zL[codes[idx]]
            px_before[idx] = new_px
            py_before[idx] = new_py
            idx = np.flatnonzero((x_before != -1) & (y_before != -1) & (sh < 0))
            new_x = x_before[idx]*cos_zR[codes[idx]] - y_before[idx]*sin_zR[codes[idx]]
            new_y = x_before[idx]*sin_zR[codes[idx]] + y_before[idx]*cos_zR[codes[idx]]
            x_before[idx] = new_x
            y_before[idx] = new_y
            idx = np.flatnonzero((px_before != -1) & (py_before != -1) & (sh < 0))
            new_px = px_before[idx]*cos_zR[codes[idx]] - py_before[idx]*sin_zR[codes[idx]]
            new_py = px_before[idx]*sin_zR[codes[idx]] + py_before[idx]*cos_zR[codes[idx]]
            px_before[idx] = new_px
            py_before[idx] = new_py
            idx = np.flatnonzero((x_after != -1) & (y_after != -1) & (sh >= 0))
            new_x = x_after[idx]*cos_zL[codes[idx]] - y_after[idx]*sin_zL[codes[idx]]
            new_y = x_after[idx]*sin_zL[codes[idx]] + y_after[idx]*cos_zL[codes[idx]]
            x_after[idx] = new_x
            y_after[idx] = new_y
            idx = np.flatnonzero((px_after != -1) & (py_after != -1) & (sh >= 0))
            new_px = px_after[idx]*cos_zL[codes[idx]] - py_after[idx]*sin_zL[codes[idx]]
            new_py = px_after[idx]*sin_zL[codes[idx]] + py_after[idx]*cos_zL[codes[idx]]
            px_after[idx] = new_px
            py_after[idx] = new_py
            idx = np.flatnonzero((x_after != -1) & (y_after != -1) & (sh < 0))
            new_x = x_after[idx]*cos_zR[codes[idx]] - y_after[idx]*sin_zR[codes[idx]]
            new_y = x_after[idx]*sin_zR[codes[idx]] + y_after[idx]*cos_zR[codes[idx]]
            x_after[idx] = new_x
            y_after[idx] = new_y
            idx = np.flatnonzero((px_after != -1) & (py_after != -1) & (sh < 0))
            new_px = px_after[idx]*cos_zR[codes[idx]] - py_after[idx]*sin_zR[codes[idx]]
            new_py = px_after[idx]*sin_zR[codes[idx]] + py_after[idx]*cos_zR[codes[idx]]
            px_after[idx] = new_px
            py_after[idx] = new_py

            # Correct for length_front as in FlukaCollimator
            idx = np.flatnonzero(s_before != -1)
            s_before[idx] = s_before[idx] - length_front[codes[idx]]
            idx = np.flatnonzero(s_after  != -1)
            s_after[idx]  = s_after[idx] - length_front[codes[idx]]

        if frame != 'jaw':
            df["s_before"]  = s_before
            df["x_before"]  = x_before
            df["px_before"] = px_before
            df["y_before"]  = y_before
            df["py_before"] = py_before
            df["s_after"]   = s_after
            df["x_after"]   = x_after
            df["px_after"]  = px_after
            df["y_after"]   = y_after
            df["py_after"]  = py_after

        return df

    # TODO: list of impacted collimators


    # TODO: does not work when multiple children
    def interactions_per_collimator(self, collimator=0, *, turn=None):
        self._check_columns_recorded(['particle_id_before'], 'interactions_per_collimator')
        n_rows = self._index.num_recorded
        if isinstance(collimator, str):
            collimator = self._collimator_id(collimator)
        mask = (self._inter[:n_rows] > 0) & (self.at_element[:n_rows] == collimator)
        if turn is not None:
            mask = mask & (self.at_turn[:n_rows] == turn)
            df = pd.DataFrame({
                    'int':  [shortcuts[inter] for inter in self._inter[:n_rows][mask]],
                    'pid':  self.particle_id_before[:n_rows][mask]
                })
            return df.groupby('pid', sort=False, group_keys=False)['int'].agg(list)
        else:
            df = pd.DataFrame({
                    'int':   [shortcuts[inter] for inter in self._inter[:n_rows][mask]],
                    'turn':  self.at_turn[:n_rows][mask],
                    'pid':   self.particle_id_before[:n_rows][mask]
                })
            return df.groupby(['pid', 'turn'], sort=False, group_keys=False)['int'].apply(list)

    def first_touch_per_turn(self, frame=None):
        self._check_columns_recorded(['particle_id_before'], 'first_touch_per_turn')
        n_rows = self._index.num_recorded
        df = pd.DataFrame({'particle_id_before': self.particle_id_before[:n_rows],
                           'at_turn': self.at_turn[:n_rows],
                           'at_element': self.at_element[:n_rows]})
        mask = np.char.startswith(self.interaction_type[:n_rows], 'Enter Jaw')
        idx_first = [group.at_element.idxmin() for _, group in df[mask].groupby(
                        ['at_turn', 'particle_id_before'], sort=False, group_keys=False)]
        df_first = self.to_pandas(frame=frame).loc[idx_first]
        df_first.insert(2, "jaw", df_first.interaction_type.astype(str).str[-1])
        to_drop = ['interaction_type',
                   *[col for col in df_first.columns if col.endswith('_after')]]
        to_rename = {col: col.replace('_before', '') for col in df_first.columns if col.endswith('before')}
        to_rename['particle_id_before'] = 'pid'
        return df_first.drop(columns=to_drop).rename(columns=to_rename)


    def _column_is_recorded(self, column):
        if column == '_index':
            return True
        if hasattr(self, '_recorded_columns'):
            return column in self._recorded_columns
        return len(getattr(self, column)) > 0

    def _check_columns_recorded(self, columns, frame):
        missing = [col for col in columns if not self._column_is_recorded(col)]
        if missing:
            raise ValueError(
                f"Cannot convert InteractionRecord to {frame} frame because "
                f"columns {missing} were not recorded.")

    def _collimator_name(self, element_id):
        if hasattr(element_id, 'tolist'):
            element_id = element_id.tolist()
        if not hasattr(self, '_coll_names'):
            return element_id
        elif element_id not in self._coll_names:
            raise ValueError(f"Element {element_id} not found in list of "
                             f"collimators of this record table!\nDid the line"
                             f" change without updating the list in the table?"
                             f"\nPlease only initialise the InteractionRecord "
                             f"after all line manipulations are complete.")
        else:
            return self._coll_names[element_id]

    def _collimator_id(self, element_name):
        if hasattr(element_name, 'tolist'):
            element_name = element_name.tolist()
        if not hasattr(self, '_coll_ids'):
            return element_name
        elif element_name not in self._coll_ids:
            raise ValueError(f"Element {element_name} not found in list of "
                             f"collimators of this record table!\nDid the line"
                             f" change without updating the list in the table?"
                             f"\nPlease only initialise the InteractionRecord "
                             f"after all line manipulations are complete.")
        else:
            return self._coll_ids[element_name]


def _make_columns(*, num_rows, columns):
    array_field_names = [ff.name for ff in InteractionRecord._XoStruct._fields
                            if hasattr(ff.ftype, 'to_nplike')]
    obligatory_columns = {'at_turn', 'at_element', 'shape_id', '_inter'}

    if columns is None:
        selected_columns = set(array_field_names)

    else:
        if not hasattr(columns, '__iter__') or isinstance(columns, str):
            columns = [columns]
        selected_columns = set(columns)
        selected_columns.discard('_index')
        selected_columns.discard('_record_all_columns')

    final_columns = []
    unknown_columns = []
    for col in selected_columns:
        if col not in array_field_names:
            if f'{col}_before' in array_field_names \
            and f'{col}_after' in array_field_names:
                final_columns.append(f'{col}_before')
                final_columns.append(f'{col}_after')
            else:
                unknown_columns.append(col)
        else:
            final_columns.append(col)
    final_columns = set(final_columns)
    unknown_columns = set(unknown_columns)

    if unknown_columns:
        raise ValueError(f"Unknown InteractionRecord columns: "
                            f"{sorted(unknown_columns)}")

    final_columns |= obligatory_columns
    init_dict = {
        field: num_rows if field in final_columns else 0
        for field in array_field_names
    }
    _record_all_columns = int(final_columns == set(array_field_names))
    _recorded_columns = tuple([field for field in array_field_names
                               if field in final_columns])

    return init_dict, _record_all_columns, _recorded_columns


def _set_recording_flags(elements, record_impacts, record_exits,
                         record_scatterings):
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
        if names is None:
            names = []
            for ii, ee in enumerate(elements):
                if hasattr(ee, 'name'):
                    names.append(ee.name)
                else:
                    name = f"el_{ii}"
                    ee.name = name
                    names.append(name)
        else:
            for nn, ee in zip(names, elements):
                if hasattr(ee, 'name'):
                    if ee.name != nn:
                        print(f"Warning: Element {nn} has name {ee.name}, but "
                              f"InteractionRecord knows this element as {nn}!")
    else:
        if elements is not None:
            raise ValueError("Cannot provide both line and elements!")
        if names is None or names is True:
            tt = line.get_table()
            names = []
            for cc in block_classes:
                ttcc = tt.rows.match(element_type=cc.__name__)
                names += list(ttcc.name)
            elements = [line.get(nn) for nn in names]
            if len(names) == 0:
                raise ValueError("No Xcoll elements in line!")
        elif names is False:
            names = []
            elements = []
        else:
            for name in names:
                if name not in line.element_names:
                    raise ValueError(f"Element {name} not found in line!")
            elements = [line.get(nn) for nn in names]
    for nn, ee in zip(names, elements):
        if not isinstance(ee, block_classes):
            raise ValueError(f"Element {nn} not an Xcoll element (expected one"
                             f" of {block_classes}, got {type(ee)})!")
    return elements, names
