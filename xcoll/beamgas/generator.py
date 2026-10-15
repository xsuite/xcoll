# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2026.                 #
# ######################################### #

import numpy as np

import xtrack as xt

from ..beam_elements.beamgas import BeamGasMarker

# Tolerance on the position of the last marker, at the end of the line [m]
_S_TOL = 1e-6


def _is_sequence(value):
    return isinstance(value, (list, tuple, np.ndarray))


def _per_marker(value, n_markers, name):
    """
    One value per marker: a sequence is checked, a scalar is repeated.
    """
    if _is_sequence(value):
        value = list(value)
        if len(value) != n_markers:
            raise ValueError(
                f"`{name}` has {len(value)} entries, but {n_markers} "
                "markers are installed.")
        return value
    return [value] * n_markers


def _is_free(env, name):
    """
    Whether `name` can be given to a new element of the environment.
    """
    if name in env.elements:
        return False
    try:
        env._check_name_clashes(name)
    except ValueError:
        return False
    return True


def _free_names(env, prefix, n_names):
    """
    The `n_names` names ``f'{prefix}.{i}'`` with the lowest free indices.
    """
    names = []
    ii = 0
    while len(names) < n_names:
        name = f'{prefix}.{ii}'
        if _is_free(env, name):
            names.append(name)
        ii += 1
    return names


class BeamGasGenerator:
    """
    Beam-gas scattering simulation of a line, based on the
    :class:`xcoll.BeamGasMarker` elements installed in it.

    An instance bound to a line is available as ``line.xcoll.beamgas``. The
    generator holds no state other than the line: the configuration is stored
    in the markers, so that a line saved to JSON after it has been configured
    generates the same particles once loaded again.

    Parameters
    ----------
    line : xtrack.Line
        Line containing, or receiving, the beam-gas markers.
    """

    def __init__(self, line):
        if not isinstance(line, xt.Line):
            raise TypeError(
                f"`line` must be an xtrack.Line (got {type(line).__name__}).")
        self._line = line

    @property
    def line(self):
        """
        Line of the beam-gas simulation.

        Returns
        -------
        line : xtrack.Line
            Line bound to this generator.
        """
        return self._line

    @property
    def marker_names(self):
        """
        Names of the beam-gas markers of the line.

        Returns
        -------
        names : list of str
            Names of the :class:`xcoll.BeamGasMarker` elements, in line order.
        """
        element_dict = self.line.element_dict
        names = [nn for nn in self.line.element_names
                 if isinstance(element_dict[nn], BeamGasMarker)]
        if len(set(names)) != len(names):
            repeated = sorted({nn for nn in names if names.count(nn) > 1})
            raise ValueError(
                f"Beam-gas markers {repeated} appear more than once in the "
                "line: each marker must be a distinct element.")
        return names

    def install_markers(self, at, *, from_=None, anchor=None,
                        from_anchor=None, names=None, prefix='beamgas',
                        close_ring=True):
        """
        Install beam-gas markers in the line.

        Each marker is a new :class:`xcoll.BeamGasMarker` element, placed with
        ``env.place`` and inserted with ``line.insert``: ``at``, ``from_``,
        ``anchor`` and ``from_anchor`` have the same meaning as there.
        Installing the markers by hand in the same way is equivalent.

        Each marker represents the section of the line from the previous
        marker (from the start of the line for the first one) to itself, so
        markers are best placed densely where the gas density or the optics
        change.

        Parameters
        ----------
        at : float, str, xdeps reference, or sequence of them
            Position of each marker, evaluated at installation: an ``s``
            position [m], an expression of the environment variables, the
            name of an element (to place the marker at its centre), or a
            reference ``'name@anchor'`` with ``anchor`` one of ``'start'``,
            ``'center'`` and ``'end'``. With ``from_``, the offset from the
            reference element.
        from_ : str or sequence of str, optional
            Element from which ``at`` is measured. ``'name@anchor'`` selects
            the anchor of the element.
        anchor : {'start', 'center', 'end'} or sequence of them, optional
            Anchor of the marker placed at ``at``. Markers are thin, so this
            only matters to follow the conventions of ``env.place``.
        from_anchor : {'start', 'center', 'end'} or sequence of them, optional
            Anchor of ``from_`` from which ``at`` is measured.
        names : str or sequence of str, optional
            Names of the new markers. By default ``f'{prefix}.{i}'``, with the
            lowest indices whose names are free in the environment.
        prefix : str, optional
            Prefix of the default names.
        close_ring : bool, optional
            If ``True`` (default), a marker is added at the end of the line
            when, after the insertion, the last marker is not there: otherwise
            the part of the line after it would belong to no section and
            :meth:`configure` refuses the line.

        Returns
        -------
        names : list of str
            Names of the installed markers, the closing one included.
        """
        line = self.line
        env = line.env

        at = list(at) if _is_sequence(at) else [at]
        n_markers = len(at)
        if n_markers == 0:
            raise ValueError("`at` must give at least one position.")
        from_ = _per_marker(from_, n_markers, 'from_')
        anchor = _per_marker(anchor, n_markers, 'anchor')
        from_anchor = _per_marker(from_anchor, n_markers, 'from_anchor')

        if names is None:
            names = _free_names(env, prefix, n_markers)
        else:
            names = _per_marker(
                [names] if isinstance(names, str) else names,
                n_markers, 'names')
            if len(set(names)) != len(names):
                raise ValueError("`names` contains repeated names.")
            for nn in names:
                if not _is_free(env, nn):
                    raise ValueError(
                        f"The name '{nn}' is already used in the environment.")

        places = []
        for nn, aa, ff, an, fa in zip(names, at, from_, anchor, from_anchor):
            env.elements[nn] = BeamGasMarker()
            places.append(env.place(
                nn, at=aa, from_=ff, anchor=an, from_anchor=fa))
        line.insert(places)

        installed = list(names)
        if close_ring and not self._last_marker_at_end():
            end_name, = _free_names(env, prefix, 1)
            env.elements[end_name] = BeamGasMarker()
            line.insert(end_name, at=line.get_length())
            installed.append(end_name)

        return installed

    def _last_marker_at_end(self):
        """
        Whether the last marker sits at the end of the line.

        Returns
        -------
        at_end : bool
            ``True`` if the line has markers and the last one is at the end
            of the line.
        """
        names = self.marker_names
        if len(names) == 0:
            return False
        tab = self.line.get_table()
        return abs(float(tab['s', names[-1]]) - float(tab.s[-1])) <= _S_TOL

    def _markers(self, at=None):
        """
        Validated names of the selected markers, in line order.

        Parameters
        ----------
        at : str, sequence of str, or None
            Selected markers. ``None`` selects all the markers of the line.

        Returns
        -------
        names : list of str
            Names of the selected markers, in line order.
        """
        names = self.marker_names
        if at is None:
            if len(names) == 0:
                raise ValueError(
                    "The line does not contain any BeamGasMarker. Install "
                    "them with `line.xcoll.beamgas.install_markers(...)`.")
            return names

        if isinstance(at, str):
            selected = [at]
        elif hasattr(at, '__iter__'):
            selected = list(at)
        else:
            raise TypeError(
                "Beam-gas markers are selected by name "
                f"(got {type(at).__name__}).")
        if len(selected) == 0:
            raise ValueError("No beam-gas marker selected.")
        position = {nn: ii for ii, nn in enumerate(names)}
        for nn in selected:
            if not isinstance(nn, str):
                raise TypeError(
                    "Beam-gas markers are selected by name "
                    f"(got {type(nn).__name__}).")
            if nn in position:
                continue
            if nn in set(self.line.element_names):
                raise TypeError(
                    f"Element '{nn}' is not a BeamGasMarker "
                    f"(got {type(self.line.element_dict[nn]).__name__}).")
            raise ValueError(f"Element '{nn}' is not present in the line.")

        return sorted(set(selected), key=position.__getitem__)
