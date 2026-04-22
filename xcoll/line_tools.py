# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
from warnings import warn
from numbers import Number

import xtrack as xt
from xtrack import line

from .beam_elements import (element_classes, collimator_classes, block_classes,
                            crystal_classes)


def _iterable(obj):
    return  hasattr(obj, '__iter__') and not isinstance(obj, str)


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


class XcollLineAPI:
    def __init__(self, line):
        self._collimators = XcollCollimatorAPI(line=line)
        self._scattering = XcollScatteringAPI(line=line)

    @property
    def collimators(self):
        return self._collimators

    @property
    def scattering(self):
        return self._scattering


class XcollEnvironmentAPI:
    # TODO
    def __init__(self, environment):
        self._environment = environment


class XcollScatteringAPI(XcollLineAccessor):

    @property
    def names(self):
        # This makes sure the accessor can access the names of the
        # collimators dynamically
        return self.line.get_elements_of_type(element_classes)[1]

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
        # This makes sure the accessor can access the names of the
        # collimators dynamically
        return self.line.get_elements_of_type(collimator_classes)[1]

    @property
    def families(self):
        families = {}
        try:
            prop_families = self.family
        except AttributeError:
            return families
        else:
            for name in self.names:
                if name in prop_families:
                    if prop_families[name] not in families:
                        families[prop_families[name]] = []
                    families[prop_families[name]].append(name)
                else:
                    if 'no family' not in families:
                        families['no family'] = []
                    families['no family'].append(name)
            return families

    def __getitem__(self, name):
        # We can getitem by name or family, so we overwrite the super method
        if name in self.families:
            return XcollLineAccessor(line=self.line, names=self.families[name])
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
        if not _iterable(names) or not _iterable(elements) \
        or not _iterable(at):
            if _iterable(names):
                raise ValueError("`names` should not be a list if any of the "
                                 "other arguments is not a list.")
            if _iterable(elements):
                raise ValueError("`elements` should not be a list if any of "
                                 "the other arguments is not a list.")
            if _iterable(at):
                raise ValueError("`at` should not be a list if any of the "
                                 "other arguments is not a list.")
            names = [names]
            elements = [elements]
            at = [at for _ in range(len(names))]
            apertures = [apertures for _ in range(len(names))]
        if not _iterable(apertures):
            apertures = [apertures for _ in range(len(names))]
        names = np.array(names)
        length = np.array([coll.length for coll in elements])
        if len(length) != len(names):
            raise ValueError("Length of `elements` does not match length of "
                             "`names`.")
        if len(at) != len(names):
            raise ValueError("Length of `at` does not match length of "
                             "`names`.")
        if len(apertures) != len(names):
            raise ValueError("Length of `apertures` does not match length of "
                             "`names`.")

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
        aper_upstream   = []
        aper_downstream = []
        for s1, s2, name, aper in zip(s_start, s_end, names, apertures):
            if not need_apertures:
                aper_upstream.append(None)
                aper_downstream.append(None)
            else:
                aper1, aper2 = self.find_aperture(name, s_start=s1, s_end=s2,
                                                  aperture=aper, table=tt,
                                                  s_tol=s_tol)
                aper_upstream.append(aper1)
                aper_downstream.append(aper2)

        # Remove elements at location of collimator (by changing them
        # into markers) and add pointers to line and names to elements
        to_remove = []
        for s1, s2, name, el in zip(s_start, s_end, names, elements):
            self.prepare_space(name, s_start=s1, s_end=s2, table=tt,
                               s_tol=s_tol, to_remove=to_remove)
            el._line = self.line
            el._name = name

        # Install
        insertions = []
        to_delete = []
        env = self.line.env
        for nn, ee, ss in zip(names, elements, s_start):
            if nn in env.elements:
                # remove placeholders with the same name
                to_remove.append(nn)
                to_delete.append(nn) # only delete original elements/markers
            insertions.append(env.place(nn, at=ss, anchor='start'))

        if len(to_remove) > 0:
            # replaces it with a drift if needed
            self.line.remove(to_remove, s_tol=s_tol)

        # remove old elements/markers from environment (after placing new
        # ones to avoid issues with names)
        for nn in to_delete:
            del env.elements[nn]

        # Add new elements to environment
        for nn, ee in zip(names, elements):
            env.elements[nn] = ee

        # Apertures
        if need_apertures:
            for s1, name, aper1, aper2 in zip(s_start, names, aper_upstream,
                                              aper_downstream):
                nn1 = f'{name}_aper_upstream'
                nn2 = f'{name}_aper_downstream'
                env.elements[nn1] = aper1
                env.elements[nn2] = aper2
                insertions.append(env.place(nn1, at=name+'@start'))
                insertions.append(env.place(nn2, at=name+'@end'))
        self.line.insert(insertions, s_tol=s_tol)

    def find_aperture(self, name, *, s_start, s_end, aperture=None, table=None,
                     s_tol=1.e-6):
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
            if table is None:
                table = self.line.get_table()
            # We look for apertures, prioritising the ones closest to the
            # centre of the collimator, and from there (if absent),
            # searching towards the upstream and downstream side.
            s_centre = (s_start + s_end)/2
            # Upstream aperture
            tt = table.rows[s_start-s_tol:s_centre+s_tol:'s']
            tt_aper = tt.rows.match('Limit.*', 'element_type')
            if len(tt_aper.name) == 0:
                aper1 = None
            else:
                aper_name = tt_aper.name[-1]
                aper1 = self.line.get(aper_name)
                if aper1.transformations_active:
                    s_aper = tt.rows[aper_name].s_start[0]
                    if not np.isclose(s_aper, s_start, atol=s_tol):
                        print(f"Warning: Using aperture {aper_name} upstream "
                              f"of {name}, but transformations are present. "
                              f"Proceed with caution.")
            # Downstream aperture
            tt = table.rows[s_centre-s_tol:s_end+s_tol:'s']
            tt_aper = tt.rows.match('Limit.*', 'element_type')
            if len(tt_aper.name) == 0:
                aper2 = None
            else:
                aper_name = tt_aper.name[0]
                aper2 = self.line.get(aper_name)
                if aper2.transformations_active:
                    s_aper = tt.rows[aper_name].s_end[0]
                    if not np.isclose(s_aper, s_end, atol=s_tol):
                        print(f"Warning: Using aperture {aper_name} downstream"
                              f" of {name}, but transformations are present. "
                              f"Proceed with caution.")
            if aper1 is None and aper2 is not None:
                aper1 = aper2
            elif aper2 is None and aper1 is not None:
                aper2 = aper1
            elif aper1 is None and aper2 is None:
                raise ValueError(f"No aperture found for {name}! Please "
                                 f"provide one.")
            return aper1, aper2

    def prepare_space(self, name, *, s_start, s_end, table=None,
                      s_tol=1.e-6, to_remove=[]):
        if table is None:
            table = self.line.get_table()
        tt = table.rows[s_start-s_tol:s_end+s_tol:'s']
        for el_name, el_type in zip(tt.name[:-1], tt.element_type[:-1]):
            if el_type == 'Marker' or el_type.startswith('Drift'):
                continue
            if not el_type.startswith('Limit'):
                print(f"Warning: Removed active element {el_name} "
                    + f"at location inside collimator {name}!")
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
        crystals = self.line.get_elements_of_type(crystal_classes)[0]
        if len(crystals) == 0:
            warn("No crystals found in line to align to beam divergence.")
        for el in crystals:
            el.align_to_beam_divergence()
