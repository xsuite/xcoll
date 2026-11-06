# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
from warnings import warn
from numbers import Number
from collections import Counter

import xtrack as xt

from .beam_elements import (element_classes, collimator_classes, block_classes,
                            crystal_classes)


def _iterable(obj):
    return  hasattr(obj, '__iter__') and not isinstance(obj, str)


class XcollLineAPI:
    def __init__(self, line):
        self._line = line
        self._scattering = XcollScatteringAPI(line=line)
        self._collimators = XcollCollimatorAPI(line=line)
        # Backwards compatibility
        if not hasattr(self.line, '_scattering') or self.line._scattering is None:
            self.line._scattering = self._scattering
        if not hasattr(self.line, '_collimators') or self.line._collimators is None:
            self.line._collimators = self._collimators

    @property
    def line(self):
        return self._line

    @property
    def scattering(self):
        """
        Interface to Xcoll scattering tools for this line.

        Returns
        -------
        scattering : object
            Xcoll scattering API bound to this line.
        """
        return self._scattering

    @property
    def collimators(self):
        """
        Interface to Xcoll collimator tools for this line.

        Returns
        -------
        collimators : object
            Xcoll collimator API bound to this line.
        """
        return self._collimators


class XcollLineAccessor:
    _typename = 'element'

    def __init__(self, line, names=None):
        super().__setattr__('_line', line)
        super().__setattr__('_env', line.env)
        if names:
            super().__setattr__('names', names)
            # self.names = names

    @property
    def env(self):
        return self._env

    @property
    def line(self):
        return self._line

    @property
    def _coll_dict(self):
        return {name: self.line.get(name) for name in self.names}

    def keys(self):
        return self._coll_dict.keys()

    def values(self):
        return self._coll_dict.values()

    def items(self):
        return self._coll_dict.items()

    def __iter__(self):
        super().__setattr__('_iter_names', iter(self.names))
        return self

    def __next__(self):
        try:
            name = next(self._iter_names)
        except StopIteration:
            raise StopIteration
        else:
            return self.line[name]

    def __len__(self):
        return len(self.names)

    def __contains__(self, key):
        return key in self.names

    def __getattr__(self, attr):
        properties = {}
        for name, el in self.items():
            if hasattr(el, attr):
                properties[name] = getattr(el, attr)
        if len(properties) == 0:
            raise AttributeError(f"Attribute `{attr}` not found.")
        if len({tuple(ii) if isinstance(ii, list) else ii
                for ii in properties.values()}) == 1:
            # If all values are the same, return a single value
            return next(iter(properties.values()))
        return properties

    def __setattr__(self, attr, value):
        if isinstance(value, dict):
            for name, el in self.items():
                if name in value:
                    if not hasattr(el, attr):
                        raise AttributeError(f"Attribute `{attr}` not found in"
                                             f" {self._typename} `{name}`.")
                    setattr(el, attr, value[name])
        else:
            # If value is not a dict, we assume it is a single value to
            # set for all collimators
            for name, el in self.items():
                if hasattr(el, attr):
                    setattr(el, attr, value)

    def __getitem__(self, name):
        if name in self.names:
            return self.line[name]
        else:
            raise ValueError(f"{self._typename.capitalize()} `{name}` not "
                             f"found in line!")

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))}>"

    def __str__(self):
        if len(self.names) == 0:
            return ''
        res = [f'{self._typename.capitalize()}s:']
        for name in self.names:
            cls_name = self.line[name].__class__.__name__
            res.append(f"    {name:<16} ({cls_name})")
        return "\n".join(res)


class XcollScatteringAPI(XcollLineAccessor):

    @property
    def names(self):
        return _element_names_of_type(self.line, element_classes)

    def enable(self):
        if len(self) == 0:
            print("No xcoll elements found in line.")
        else:
            nemitt_x = None
            nemitt_y = None
            for el in self:
                if hasattr(el, 'optics') and el.optics is not None:
                    if nemitt_x is None:
                        nemitt_x = el.nemitt_x
                    if nemitt_y is None:
                        nemitt_y = el.nemitt_y
                    if not np.isclose(el.nemitt_x, nemitt_x) \
                    or not np.isclose(el.nemitt_y, nemitt_y):
                        raise ValueError("Not all collimators have the same "
                                         "emittance. This is not supported.")
                if hasattr(el, 'enable_scattering'):
                    el.enable_scattering()

    def disable(self):
        if len(self) == 0:
            print("No xcoll elements found in line.")
        else:
            for el in self:
                if hasattr(el, 'disable_scattering'):
                    el.disable_scattering()

    def identify_primary_losses(self):
        if len(self) == 0:
            print("No xcoll elements found in line.")
        else:
            for el in self:
                el.mark_scattered_particles = True


class XcollCollimatorAPI(XcollLineAccessor):
    _typename = 'collimator'

    def __str__(self):
        res = []
        if len(self.families) > 0:
            res.append('Families:')
            for family, names in self.families.items():
                res.append(f"    {family:8}: {', '.join(names)}")
            res.append('')
        res.append(super().__str__())
        return "\n".join(res)

    @property
    def names(self):
        return _element_names_of_type(self.line, collimator_classes)

    @property
    def families(self):
        families = {}
        missing = object()
        has_family = False
        for name in self.names:
            family = getattr(self.line.get(name), 'family', missing)
            if family is missing:
                family = 'no family'
            else:
                has_family = True
            families.setdefault(family, []).append(name)
        return families if has_family else {}

    def __getitem__(self, name):
        # We can getitem by name or family, so we overwrite the super method
        families = self.families
        if name in families:
            return XcollLineAccessor(line=self.line, names=families[name])
        elif name in self.names:
            return self.line[name]
        else:
            raise ValueError(f"Neither family nor collimator `{name}` found in"
                             f" line!")

    def open(self, names=None):
        if names is None:
            names = self.names
        if len(names) == 0:
            print("No collimators found in line.")
        else:
            for coll in names:
                self.line.get(coll).open_jaws(keep_tilts=False)
                self.line.get(coll).gap = None

    def to_parking(self, names=None):
        if names is None:
            names = self.names
        if len(names) == 0:
            print("No collimators found in line.")
        else:
            raise NotImplementedError("Need to move this to new type manager "
                                      "or so.")

    def install(self, names, elements, *, at=None, apertures=None,
                need_apertures=False, s_tol=1.e-6, at_s=None):
        if at_s is not None:
            warn("Warning: `at_s` is deprecated and will be removed in the "
                 "future. Please use `at` instead.", FutureWarning)
            at = at_s
        if not _iterable(names):
            names = [names]
        if not _iterable(elements):
            elements = [elements]
        if not _iterable(at):
            at = [at for _ in range(len(names))]
        if not _iterable(apertures) or len(apertures) != len(names):
            # A single aperture spec (None, an aperture name, a single
            # aperture object, or an [upstream, downstream] pair) is applied
            # to every element. Only when `apertures` already provides
            # exactly one entry per name is it used as a genuine per-element
            # list (one aperture spec per collimator). Each of those entries
            # can in turn be a single aperture spec or its own
            # [upstream, downstream] pair, e.g. for 3 collimators:
            #     apertures = [ap1, [ap2_up, ap2_down], 'ap3_name']
            apertures = [apertures for _ in range(len(names))]
        # names = np.array(names)
        if len(elements) != len(names):
            raise ValueError("Length of `elements` does not match length of "
                             "`names`.")
        if len(at) != len(names):
            raise ValueError("Length of `at` does not match length of "
                             "`names`.")
        if len(apertures) != len(names):
            raise ValueError("Length of `apertures` does not match length of "
                             "`names`.")
        length = np.array([coll.length for coll in elements])

        # Verify elements
        for name, el in zip(names, elements):
            if not isinstance(el, block_classes):
                raise ValueError(f"Element {el} is not a valid block or "
                                 f"collimator class.")
            el._tracking = False
            if el.name is None:
                el.name = name

        # Get positions
        tt = self.line.get_table()
        s_start = []
        for name, s, l in zip(names, at, length):
            if s is not None and not isinstance(s, Number):
                # TODO: this could be generalised to allow the same API as
                # line.insert, however, this has to be implemented cautiously
                # and well-tested.
                raise NotImplementedError("`at` must be a number indicating "
                                    "the s position of the blow-up element.")
            existing_length = 0
            if name in tt.name:
                if hasattr(self.line[name], 'length'):
                    existing_length = self.line[name].length
                else:
                    existing_length = 0
                existing_s = tt.rows[name].s_start[0]
                new_s = existing_s + existing_length/2. - l/2
                if s is not None and not np.isclose(s, new_s, atol=s_tol):
                    raise ValueError(f"Element {name} already exists in line "
                            f"at location {existing_s} with length "
                            f"{existing_length}. Provided `at` = {s} is "
                            f"different from the existing one. Please provide "
                            f"an `at` for this element that corresponds to the"
                            f" existing location or use `at=None`, or remove "
                            f"the element first.")
                s = new_s
            elif s is None:
                raise ValueError(f"Element {name} not found in line. Need to "
                                 f"manually provide `at`.")
            s_start.append(s)
        s_start = np.array(s_start)
        s_end = s_start + length

        # Check positions
        l_line = tt.s_end[-1]
        for s1, s2, name in zip(s_start, s_end, names):
            if s1 <= s_tol:
                raise ValueError(f"Position of {name} too close to start of "
                                 f"line. Please cycle.")
            if s2 >= l_line - s_tol:
                raise ValueError(f"Position of {name} too close to end of "
                                 f"line. Please cycle.")

        # Look for apertures
        aper_upstream = []
        aper_downstream = []
        if need_apertures:
            tt_apertures = tt.rows.match('Limit.*', 'element_type')

        for s1, s2, name, aper in zip(s_start, s_end, names, apertures):
            if not need_apertures:
                aper_upstream.append(None)
                aper_downstream.append(None)
            else:
                aper1, aper2 = self.find_aperture(
                    name, s_start=s1, s_end=s2,
                    aperture=aper, aperture_table=tt_apertures,
                    s_tol=s_tol
                )
                aper_upstream.append(aper1)
                aper_downstream.append(aper2)

        # Identify collimators that can replace existing elements in place
        direct_replacements = set()
        occurrences = Counter(self.line.element_names)
        drift_like_classes = (xt.Drift, xt.DriftExact, xt.Device)
        for name, el, s1, s2 in zip(names, elements, s_start, s_end):
            if occurrences[name] != 1:
                continue
            old = self.line.get(name)
            if not isinstance(old, drift_like_classes):
                continue
            if not getattr(old, 'isthick', False):
                continue
            idx = tt.rows.get_index(name)
            old_length = tt.s_end[idx] - tt.s_start[idx]
            if not np.isclose(old_length, el.length, atol=s_tol, rtol=0):
                continue
            if not np.isclose(tt.s_start[idx], s1, atol=s_tol, rtol=0):
                continue
            direct_replacements.add(name)

        # Remove elements at location of collimator (by changing them
        # into markers) and add pointers to line and names to elements
        to_remove = []
        for s1, s2, name, el in zip(s_start, s_end, names, elements):
            if name not in direct_replacements:
                self.prepare_space(
                    name, s_start=s1, s_end=s2, table=tt, s_tol=s_tol,
                    to_remove=to_remove
                )
            el._line = self.line
            el._name = name

        # Install
        env = self.line.env
        insertions = []
        to_delete = []
        # Only the fallback elements require normal insertion
        for nn, ee, ss in zip(names, elements, s_start):
            if nn in direct_replacements:
                continue
            if nn in env.elements:
                to_remove.append(nn)
                to_delete.append(nn)
            insertions.append(env.place(nn, at=ss, anchor='start'))

        # We are modifying the line, even if no insert/remove is needed
        self.line._method_incompatible_with_compose()
        self.line.discard_tracker()

        if to_remove:
            # Replaces it with a drift if needed
            self.line.remove(list(dict.fromkeys(to_remove)), s_tol=s_tol)

        # Remove the old placeholders for fallback installations
        for nn in to_delete:
            del env.elements[nn]

        # Register fallback collimators
        for nn, ee in zip(names, elements):
            if nn not in direct_replacements:
                env.elements[nn] = ee

        # Insert only collimators requiring geometric changes
        if insertions:
            self.line.insert(insertions, s_tol=s_tol)

        # Directly replace existing drift-like placeholders
        for nn, ee in zip(names, elements):
            if nn in direct_replacements:
                env._element_dict[nn] = ee

        # Add apertures directly, without calling line.insert()
        if need_apertures:
            aperture_names = {}
            for name, aper1, aper2 in zip(
                names, aper_upstream, aper_downstream
            ):
                nn1 = f'{name}_aper_upstream'
                nn2 = f'{name}_aper_downstream'
                env.elements[nn1] = aper1
                env.elements[nn2] = aper2
                aperture_names[name] = (nn1, nn2)
            new_names = []
            for nn in self.line.element_names:
                if nn in aperture_names:
                    # Install apertures immediately before/after each collimator.
                    upstream, downstream = aperture_names[nn]
                    new_names.extend((upstream, nn, downstream))
                else:
                    new_names.append(nn)
            self.line.element_names[:] = new_names

    def find_aperture(self, name, *, s_start, s_end, aperture=None,
                      table=None, aperture_table=None, s_tol=1.e-6):
        if aperture is not None:
            if isinstance(aperture, str):
                try:
                    aper1 = self.env.get(aperture)
                except KeyError:
                    aper1 = self.line.get(aperture)
                aper2 = aper1
            elif _iterable(aperture):
                if len(aperture) == 1:
                    aperture = [aperture[0], aperture[0]]
                elif len(aperture) != 2:
                    raise ValueError("The value `aperture` should be None or "
                                     "a list [upstream, downstream].")
                if aperture[0] is None or aperture[1] is None:
                    raise ValueError("Both upstream and downstream apertures "
                                     "must be provided if `aperture` is a "
                                     "list.")
                if isinstance(aperture[0], str):
                    try:
                        aper1 = self.env.get(aperture[0])
                    except KeyError:
                        aper1 = self.line.get(aperture[0])
                else:
                    aper1 = aperture[0]
                if isinstance(aperture[1], str):
                    try:
                        aper2 = self.env.get(aperture[1])
                    except KeyError:
                        aper2 = self.line.get(aperture[1])
                else:
                    aper2 = aperture[1]
            else:
                aper1 = aperture
                aper2 = aperture
            if not xt.line._is_aperture(aper1, self.line):
                raise ValueError(f"Not a valid aperture: {aper1}")
            if not xt.line._is_aperture(aper2, self.line):
                raise ValueError(f"Not a valid aperture: {aper2}")
            return aper1, aper2

        else:
            if aperture_table is None:
                if table is None:
                    table = self.line.get_table()
                aperture_table = table.rows.match('Limit.*', 'element_type')
            s = aperture_table.s
            s_centre = (s_start + s_end) / 2
            # Upstream aperture: last aperture between start and centre
            i0 = np.searchsorted(s, s_start - s_tol, side='left')
            i1 = np.searchsorted(s, s_centre + s_tol, side='right')
            aper1 = None
            if i0 < i1:
                idx = i1 - 1
                aper_name = aperture_table.name[idx]
                aper1 = self.line.get(aper_name)
                if aper1.transformations_active:
                    s_aper = aperture_table.s_start[idx]
                    if not np.isclose(s_aper, s_start, atol=s_tol):
                        print(f"Warning: Using aperture {aper_name} upstream "
                            f"of {name}, but transformations are present. "
                            f"Proceed with caution.")
            # Downstream aperture: first aperture between centre and end
            i0 = np.searchsorted(s, s_centre - s_tol, side='left')
            i1 = np.searchsorted(s, s_end + s_tol, side='right')
            aper2 = None
            if i0 < i1:
                idx = i0
                aper_name = aperture_table.name[idx]
                aper2 = self.line.get(aper_name)
                if aper2.transformations_active:
                    s_aper = aperture_table.s_end[idx]
                    if not np.isclose(s_aper, s_end, atol=s_tol):
                        print(f"Warning: Using aperture {aper_name} downstream "
                            f"of {name}, but transformations are present. "
                            f"Proceed with caution.")
            # Fallback if only one aperture was found
            if aper1 is None and aper2 is not None:
                aper1 = aper2
            elif aper2 is None and aper1 is not None:
                aper2 = aper1
            elif aper1 is None and aper2 is None:
                raise ValueError(f"No aperture found for {name}! "
                                f"Please provide one.")
            return aper1, aper2

    def prepare_space(self, name, *, s_start, s_end, table=None,
                      s_tol=1.e-6, to_remove=None):
        if table is None:
            table = self.line.get_table()
        if to_remove is None:
            to_remove = []
        # Find the relevant rows without constructing a new table
        s = table.s
        i0 = np.searchsorted(s, s_start - s_tol, side='left')
        i1 = np.searchsorted(s, s_end + s_tol, side='right')
        # Preserve the original exclusion of the final selected row
        i1 = max(i0, i1 - 1)
        for el_name, el_type in zip(table.name[i0:i1], table.element_type[i0:i1]):
            if (el_type in ('Marker', 'Device', 'ThickSliceDevice')
                or el_type.startswith('Drift')):
                continue
            if not el_type.startswith('Limit'):
                print(f"Warning: Removed active element {el_name} "
                    f"at location inside collimator {name}!")
            to_remove.append(el_name)
        return to_remove

    def get_optics_at(self, names, *, twiss=None):
        if twiss is None:
            twiss = self.line.twiss()
        if not _iterable(names):
            names = [names]
        coll_entry_indices = twiss.rows.indices[names]
        tw_entry = twiss.rows[coll_entry_indices]
        tw_exit = twiss.rows[coll_entry_indices+1]
        tw_exit.name = tw_entry.name
        return tw_entry, tw_exit

    def assign_optics(self, *, nemitt_x=None, nemitt_y=None, twiss=None):
        tw_upstream, tw_downstream = self.get_optics_at(self.names, twiss=twiss)
        beta_gamma_rel = self.line.particle_ref._xobject.gamma0[0]*self.line.particle_ref._xobject.beta0[0]
        for name, coll in self.items():
            coll.assign_optics(name=name, nemitt_x=nemitt_x, nemitt_y=nemitt_y, twiss_upstream=tw_upstream,
                               twiss_downstream=tw_downstream, beta_gamma_rel=beta_gamma_rel)

    def align_to_beam_divergence(self):
        tt = self.line.get_table()
        crystal_names = []
        for cc in crystal_classes:
            ttcc = tt.rows.match(element_type=cc.__name__)
            crystal_names += list(ttcc.name)
        crystals = [self.line.get(nn) for nn in crystal_names]
        if len(crystals) == 0:
            warn("No crystals found in line to align to beam divergence.")
        for el in crystals:
            el.align_to_beam_divergence()


def _element_names_of_type(line, classes):
    # Quick accessor to not pass via line.get_table() which is slow
    # Preserve the existing grouping by element class
    groups = {cls.__name__: [] for cls in classes}
    elements = line.element_dict

    for name in line.element_names:
        el = elements[name]
        if isinstance(el, xt.Replica):
            el = el.resolve(line)

        group = groups.get(type(el).__name__)
        if group is not None:
            group.append(name)

    return [name for group in groups.values() for name in group]
