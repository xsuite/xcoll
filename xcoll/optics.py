# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

from collections.abc import Mapping
from numbers import Real

import numpy as np
import xtrack as xt
import xdeps as xd


_DEFAULT_POSITIONS = ('upstream', 'centre', 'downstream')


def _position_offset(position, lengths):
    if position == 'upstream':
        return np.zeros_like(lengths)
    if position == 'centre':
        return lengths / 2
    if position == 'downstream':
        return lengths
    if isinstance(position, Real) and not isinstance(position, (bool, np.bool_)):
        if not np.isfinite(position) or position < 0 or np.any(position > lengths):
            raise ValueError(f"Position {position} is outside the collimator.")
        return np.full_like(lengths, position)
    raise ValueError(f"Invalid optics position: {position!r}.")


class CollimatorOpticsTable(xt.Table):
    """Reference optics, indexed by collimator name or local position."""

    def copy(self):
        return self.rows[np.arange(len(self))]

    def __getitem__(self, key):
        if isinstance(key, str) and key in ('upstream', 'downstream') and self._index == 'position':
            from warnings import warn
            warn(f"optics[{key!r}] is deprecated; use reference_optics.rows[{key!r}] "
                 "or reference_optics[field, position] instead.", FutureWarning, stacklevel=2)
            return self.rows[key]
        return super().__getitem__(key)


_FIELDS = ('betx', 'bety', 'alfx', 'alfy', 'gamx', 'gamy', 'mux', 'muy',
           'x', 'px', 'y', 'py', 'zeta', 'delta', 'dx', 'dpx', 'dy', 'dpy')
_METADATA = ('beta_gamma_rel', 'nemitt_x', 'nemitt_y')


def _get_optics_at_positions(twiss, names, lengths, positions):
    """One table for all names; s is absolute and s_rel is local."""
    names = list(names)
    lengths = np.asarray(lengths, dtype=float)
    if isinstance(positions, (str, Real)):
        positions = [positions]
    positions = list(positions)
    if not positions or len(set(positions)) != len(positions):
        raise ValueError("Provide nonempty, distinct optics positions.")
    if len(set(names)) != len(names):
        raise ValueError("Provide distinct collimator names.")
    offsets = [_position_offset(position, lengths) for position in positions]
    indices = twiss.rows.indices[names]
    up = twiss.rows[indices]
    down = twiss.rows[indices + 1]
    if not np.allclose(down.s - up.s, lengths, rtol=1e-7, atol=1e-9):
        raise ValueError("Downstream twiss not compatible with collimator lengths.")
    rows = [_drift_optics(up, down, offset, lengths) for offset in offsets]
    for position, row in zip(positions, rows):
        if position == 'downstream':
            for field in row._col_names:
                if field not in ('name', 's') and field in down._col_names:
                    row[field] = down[field]
    # Flatten collimator-major, sorted by local s within each collimator.
    rel = np.stack(offsets, axis=1)
    order = np.argsort(rel, axis=1, kind='stable')
    labels = np.array([p if isinstance(p, str) else f's={float(p):.16g}'
                       for p in positions], dtype=object)
    def flatten(a):
        return np.take_along_axis(a, order, axis=1).ravel()
    data = dict(name=np.repeat(names, len(positions)),
                position=flatten(np.broadcast_to(labels, rel.shape)),
                s=flatten(up.s[:, None] + rel), s_rel=flatten(rel))
    for field in _FIELDS:
        if all(field in row._col_names for row in rows):
            data[field] = flatten(np.stack([row[field] for row in rows], axis=1))
    return CollimatorOpticsTable(data)


def _drift_optics(up, down, offset, lengths):
    # Nominal paraxial linear drift model. Use gamma from Xtrack when available (the
    # uncoupled identity is not valid for general coupled normal modes).
    # Do not copy unrelated Twiss columns (e.g. W_matrix or phase advances)
    # to an interior location: those would not be valid there.
    fraction = np.divide(offset, lengths, out=np.zeros_like(offset), where=lengths != 0)
    data = {'name': up.name.copy(), 's': up.s + offset}
    for plane in ('x', 'y'):
        beta = up[f'bet{plane}']
        alpha = up[f'alf{plane}']
        gamma = (up[f'gam{plane}'] if f'gam{plane}' in up._col_names
                 else (1 + alpha**2) / beta)
        at_exit = (offset == lengths) & (lengths != 0)
        data[f'bet{plane}'] = np.where(
            at_exit, down[f'bet{plane}'], beta - 2*offset*alpha + offset**2*gamma)
        data[f'alf{plane}'] = np.where(at_exit, down[f'alf{plane}'], alpha - offset*gamma)
        gamma_down = (down[f'gam{plane}'] if f'gam{plane}' in down._col_names
                      else (1 + down[f'alf{plane}']**2) / down[f'bet{plane}'])
        data[f'gam{plane}'] = np.where(at_exit, gamma_down, gamma)
    # The orbit and dispersion are affine in s in a drift. Interpolate the
    # boundary values, rather than treating canonical momenta as slopes.
    for field in ('x', 'y', 'dx', 'dy', 'px', 'py', 'dpx', 'dpy'):
        if field in up._col_names and field in down._col_names:
            data[field] = up[field] + fraction*(down[field] - up[field])
    return xt.Table(data)


def _select_optics_rows(twiss, name, length):
    """Select independent local rows; never retain the whole line/table."""
    if isinstance(twiss, Mapping):
        rows = []
        for position in ('upstream', 'downstream'):
            table = twiss[position]
            if name is None and len(table) != 1:
                raise ValueError("Need a name or single-row boundary tables.")
            row = table.rows[[0]] if name is None else table.rows[name]
            if len(row) != 1:
                raise ValueError("Need exactly one optics row for each position.")
            rows.append(row)
        start = rows[0].s[0]
        if not np.isclose(rows[1].s[0] - start, length, rtol=1e-7, atol=1e-9):
            raise ValueError("Downstream twiss not compatible with collimator length.")
        data = dict(name=np.array([name or rows[0].name[0]] * 2),
                    position=np.array(['upstream', 'downstream']),
                    s=np.array([start, start + length]), s_rel=np.array([0., length]))
        for field in _FIELDS:
            if all(field in row._col_names for row in rows):
                data[field] = np.array([row[field][0] for row in rows])
        twiss = CollimatorOpticsTable(data)
    if not isinstance(twiss, xd.Table):
        raise TypeError("Need an Xsuite optics table.")
    if 'position' not in twiss._col_names or 's_rel' not in twiss._col_names:
        if name is None:
            raise ValueError("When using a full `twiss` table, need to provide the name as well.")
        twiss = _get_optics_at_positions(twiss, [name], [length], _DEFAULT_POSITIONS)
    names = np.unique(twiss.name)
    if name is None:
        if len(names) != 1:
            raise ValueError("Need a collimator name for a combined optics table.")
        name = names[0]
    local = twiss.rows[twiss.name == name]
    if not {'upstream', 'downstream'} <= set(local.position):
        raise ValueError("Optics assignment requires upstream and downstream positions.")
    if len(set(local.position)) != len(local.position):
        raise ValueError("Local optics position labels must be unique.")
    cols = list(local._col_names)
    required = {'name', 'position', 's', 's_rel', 'betx', 'bety', 'alfx', 'alfy', 'x', 'y'}
    if not required <= set(cols):
        raise ValueError(f"Missing local optics columns: {required - set(cols)}")
    data = {key: np.array(local[key], copy=True) for key in cols}
    data.update({key: local[key] for key in local.keys(exclude_columns=True)})
    out = CollimatorOpticsTable(data, col_names=cols, index='position')
    if (not np.all(np.isfinite(out.s_rel)) or np.any(out.s_rel < -1e-9)
            or np.any(out.s_rel > length + 1e-9)
            or not np.isclose(out['s_rel', 'upstream'], 0., atol=1e-9)
            or not np.isclose(out['s_rel', 'downstream'], length, rtol=1e-7, atol=1e-9)
            or not np.allclose(out.s - out['s', 'upstream'], out.s_rel, atol=1e-9)):
        raise ValueError("Local twiss coordinates not compatible with collimator length.")
    # Per-collimator metadata may be scalar or columns on an aggregate table.
    for key in _METADATA:
        if key in out._col_names:
            values = out[key]
            if not np.all(values == values[0]):
                raise ValueError(f"Inconsistent {key} within collimator {name}.")
            data[key] = values[0]
            cols.remove(key)
    data.update(collimator_name=name, s_start=float(out['s', 'upstream']))
    return CollimatorOpticsTable(data, col_names=cols, index='position').rows[np.argsort(out.s_rel, kind='stable')]


def _combine_reference_optics(named_tables):
    """Aggregate assigned snapshots, without retaining line or element references."""
    if not named_tables:
        return CollimatorOpticsTable(dict(name=np.array([], dtype=object),
            position=np.array([], dtype=object), s=np.array([]), s_rel=np.array([])))
    _, tables = zip(*named_tables)
    cols = list(dict.fromkeys(key for table in tables for key in table._col_names))
    data = {}
    for key in cols:
        example = next(table[key] for table in tables if key in table._col_names)
        # Optional columns may differ between externally supplied optics tables.
        data[key] = np.concatenate([table[key] if key in table._col_names else
            np.full((len(table), *example.shape[1:]), np.nan) for table in tables])
    data['name'] = np.concatenate([np.full(len(table), name, dtype=object)
                                   for name, table in named_tables])
    for key in _METADATA:
        values = [table[key] for table in tables]
        if all(value == values[0] for value in values):
            data[key] = values[0]
        else:
            data[key] = np.concatenate([np.full(len(table), value) for table, value in zip(tables, values)])
            if key not in cols:
                cols.append(key)
    return CollimatorOpticsTable(data, col_names=cols)
