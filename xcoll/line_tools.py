# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
from warnings import warn

import xtrack as xt

from .beam_elements import element_classes, _all_collimator_classes, _all_block_classes


class XcollScatteringAPI:
    def __init__(self, line):
        self._line = line

    @property
    def line(self):
        return self._line

    def enable(self):
        elements = self.line.get_elements_of_type(element_classes)[0]
        if len(elements) == 0:
            print("No xcoll elements found in line.")
        else:
            nemitt_x = None
            nemitt_y = None
            for el in elements:
                if hasattr(el, 'optics') and el.optics is not None:
                    if nemitt_x is None:
                        nemitt_x = el.nemitt_x
                    if nemitt_y is None:
                        nemitt_y = el.nemitt_y
                    if not np.isclose(el.nemitt_x, nemitt_x) \
                    or not np.isclose(el.nemitt_x, nemitt_x):
                        raise ValueError("Not all collimators have the same "
                                    + "emittance. This is not supported.")
                if hasattr(el, 'enable_scattering'):
                    el.enable_scattering()

    def disable(self):
        elements = self.line.get_elements_of_type(element_classes)[0]
        if len(elements) == 0:
            print("No xcoll elements found in line.")
        else:
            for el in elements:
                if hasattr(el, 'disable_scattering'):
                    el.disable_scattering()


class XcollCollimatorAPI:
    def __init__(self, line):
        self._line = line

    @property
    def line(self):
        return self._line

    def install(self, names, elements, *, at_s=None, apertures=None, need_apertures=False, s_tol=1.e-6):
        if self.line._has_valid_tracker():
            raise Exception("Tracker already built!\nPlease install collimators before building "
                        + "tracker!")

        if not hasattr(names, '__iter__') or isinstance(names, str):
            names = [names]
        if not hasattr(elements, '__iter__') or isinstance(elements, str):
            elements = [elements]
        names = np.array(names)
        length = np.array([coll.length for coll in elements])
        assert len(length) == len(names)
        if not hasattr(at_s, '__iter__'):
            at_s = [at_s for _ in range(len(names))]
        assert len(at_s) == len(names)
        if isinstance(apertures, str) or not hasattr(apertures, '__iter__'):
            apertures = [apertures for _ in range(len(names))]
        assert len(apertures) == len(names)

        # Verify elements
        for el in elements:
            print(el.__class__)
            assert isinstance(el, _all_block_classes)
            el._tracking = False

        # Get positions
        tab = self.line.get_table()
        tt = tab.rows[[name for name in names if name in self.line.element_names]]
        s_start = []
        for name, s, l in zip(names, at_s, length):
            if s is None:
                s_start.append(self._get_s_start(name, length=l, table=tt))
            else:
                s_start.append(s)
        s_start = np.array(s_start)
        s_end = s_start + length

        # Check positions
        l_line = self.line.get_length()
        for s1, s2, name, s3 in zip(s_start, s_end, names, at_s):
            self.check_position(name, s_start=s1, s_end=s2, at_s=s3, length=l_line, s_tol=s_tol)

        # Look for apertures
        aper_upstream   = []
        aper_downstream = []
        for s1, s2, name, aper in zip(s_start, s_end, names, apertures):
            if not need_apertures:
                aper_upstream.append(None)
                aper_downstream.append(None)
            else:
                aper1, aper2 = self.get_aperture(name, s_start=s1, s_end=s2, aperture=aper, table=tab, s_tol=s_tol)
                aper_upstream.append(aper1)
                aper_downstream.append(aper2)

        # Remove elements at location of collimator (by changing them into markers)
        for s1, s2, name, el in zip(s_start, s_end, names, elements):
            self.prepare_space(name, s_start=s1, s_end=s2, table=tab, s_tol=s_tol)
            el._line = self.line
            el._name = name

        # Install
        self.line._insert_thick_elements_at_s(element_names=list(names), elements=elements, at_s=s_start, s_tol=s_tol)

        # Install apertures
        if need_apertures:
            for s1, name, aper1, aper2 in zip(s_start, names, aper_upstream, aper_downstream):
                self.line.insert_element(element=aper1, name=f'{name}_aper_upstream', at=name, s_tol=s_tol)
                idx = self.line.element_names.index(name) + 1
                self.line.insert_element(element=aper2, name=f'{name}_aper_downstream', at=idx, s_tol=s_tol)

    def check_position(self, name, *, s_start, s_end, at_s, length=None, s_tol=1.e-6):
        if at_s is None:
            if name not in self.line.element_names:
                raise ValueError(f"Element {name} not found in line. Provide `at_s`.")
        elif name in self.line.element_names:
            if at_s < s_start or at_s > s_end:
                raise ValueError(f"Element {name} already exists in line at different "
                            + f"location: at_s = {at_s}, exists at [{s_start}, {s_end}].")
        if length is None:
            length = self.line.get_length()
        if s_start <= s_tol:
            raise ValueError(f"Position of {name} too close to start of line. Please cycle.")
        if s_end >= length - s_tol:
            raise ValueError(f"Position of {name} too close to end of line. Please cycle.")

    def get_apertures_at_s(self, s, *, table=None, s_tol=1.e-6):
        if table is None:
            table = self.line.get_table()
        tab_s = table.rows[s-s_tol:s+s_tol:'s']
        aper = tab_s.rows[[cls.startswith('Limit') for cls in tab_s.element_type]]
        if len(aper) == 0:
            return None
        elif len(aper) == 1:
            return aper.name[0]
        else:
            raise ValueError(f"Multiple apertures found at location {s} with "
                           + f"tolerance {s_tol}: {aper.name}. Not supported.")

    def get_aperture(self, name, *, s_start, s_end, aperture=None, table=None, s_tol=1.e-6):
        if aperture is not None:
            if isinstance(aperture, str):
                aper1 = self.line[aperture]
                aper2 = self.line[aperture]
            elif hasattr(aperture, '__iter__'):
                if len(aperture) != 2:
                    raise ValueError(f"The value `aperture` should be None or a list "
                                + f"[upstream, downstream].")
                assert aperture[0] is not None and aperture[1] is not None
                if isinstance(aperture[0], str):
                    aper1 = self.line[aperture[0]]
                if isinstance(aperture[1], str):
                    aper2 = self.line[aperture[1]]
            else:
                aper1 = aperture
                aper2 = aperture
            if not xt.line._is_aperture(aper1, self.line):
                raise ValueError(f"Not a valid aperture: {aper1}")
            if not xt.line._is_aperture(aper2, self.line):
                raise ValueError(f"Not a valid aperture: {aper2}")
            return aper1.copy(), aper2.copy()
        else:
            if table is None:
                table = self.line.get_table()
            aper1 = self.get_apertures_at_s(s=s_start, table=table, s_tol=s_tol)
            aper2 = self.get_apertures_at_s(s=s_end, table=table, s_tol=s_tol)
            if aper1 is None and aper2 is not None:
                aper1 = aper2
                print(f"Warning: Could not find upstream aperture for {name}! "
                    + f"Used copy of downstream aperture. Proceed with caution.")
            elif aper2 is None and aper1 is not None:
                aper2 = aper1
                print(f"Warning: Could not find downstream aperture for {name}! "
                    + f"Used copy of upstream aperture. Proceed with caution.")
            elif aper1 is None and aper2 is None:
                aper_mid = self.get_apertures_at_s(s=(s_start+s_end)/2, table=table, s_tol=s_tol)
                if aper_mid is None:
                    raise ValueError(f"No aperture found for {name}! Please provide one.")
                if self.line[aper_mid].allow_rot_and_shift \
                and xt.base_element._tranformations_active(self.line[aper_mid]):
                    print(f"Warning: Using the centre aperture for {name}, but "
                        + f"transformations are present. Proceed with caution.")
                aper1 = aper_mid
                aper2 = aper_mid
            return self.line[aper1].copy(), self.line[aper2].copy()

    def prepare_space(self, name, *, s_start, s_end, table=None, s_tol=1.e-6):
        if table is None:
            table = self.line.get_table()
        tt = table.rows[s_start-s_tol:s_end+s_tol:'s']
        for element_name, element_type in zip(tt.name[:-1], tt.element_type[:-1]):
            if element_type == 'Marker' or element_type.startswith('Drift'):
                continue
            if not element_type.startswith('Limit'):
                print(f"Warning: Removed active element {element_name} "
                    + f"at location inside collimator {name}!")
            length = self.line[element_name].length if hasattr(self.line[element_name], 'length') else 0
            self.line.element_dict[element_name] = xt.Drift(length=length)

    def get_optics_at(self, names, *, twiss=None, tw=None):
        if tw is not None:
            warn("The argument tw is deprecated. Please use twiss instead.", FutureWarning)
            if twiss is None:
                twiss = tw
        if twiss is None:
            if not self.line._has_valid_tracker():
                raise Exception("Please build the tracker before computing the optics for the openings!")
            twiss = self.line.twiss()
        if not hasattr(names, '__iter__') and not isinstance(names, str):
            names = [names]
        coll_entry_indices = twiss.rows.indices[names]
        tw_entry = twiss.rows[coll_entry_indices]
        tw_exit = twiss.rows[coll_entry_indices+1]
        tw_exit.name = tw_entry.name
        return tw_entry, tw_exit

    def assign_optics(self, *, nemitt_x=None, nemitt_y=None, twiss=None, tw=None):
        if tw is not None:
            warn("The argument tw is deprecated. Please use twiss instead.", FutureWarning)
            if twiss is None:
                twiss = tw
        if not self.line._has_valid_tracker():
            raise Exception("Please build tracker before setting the openings!")
        names = self.line.get_elements_of_type(_all_collimator_classes, _all_block_classes)[1]
        tw_upstream, tw_downstream = self.get_optics_at(names, twiss=twiss)
        beta_gamma_rel = self.line.particle_ref._xobject.gamma0[0]*self.line.particle_ref._xobject.beta0[0]
        for coll in names:
            self.line[coll].assign_optics(name=coll, nemitt_x=nemitt_x, nemitt_y=nemitt_y, twiss_upstream=tw_upstream,
                                    twiss_downstream=tw_downstream, beta_gamma_rel=beta_gamma_rel)

    def open(self, names=None):
        if names is None:
            names = self.line.get_elements_of_type(_all_collimator_classes, _all_block_classes)[1]
        if len(names) == 0:
            print("No collimators found in line.")
        else:
            for coll in names:
                self.line[coll].open_jaws(keep_tilts=False)
                self.line[coll].gap = None

    def to_parking(self, names=None):
        if names is None:
            names = self.line.get_elements_of_type(_all_collimator_classes)[1]
        if len(names) == 0:
            print("No collimators found in line.")
        else:
            raise NotImplementedError("Need to move this to new type manager or so.")

    def _get_s_start(self, name, *, length, table=None):
        if table is None:
            table = self.line.get_table()
        if name in self.line.element_names and hasattr(self.line[name], 'length'):
            existing_length = self.line[name].length
        else:
            existing_length = 0
        if name not in table.name:
            raise ValueError(f"Element {name} not found in line. Need to manually provide `at_s`.")
        return table.rows[name].s[0] + existing_length/2. - length/2



# Deprecated; to be removed
# -------------------------

def assign_optics_to_collimators(line, nemitt_x=None, nemitt_y=None, twiss=None):
    warn("The function xcoll.assign_optics_to_collimators() is deprecated and will be "
       + "removed in the future. Please use line.collimators.assign_optics() instead.", FutureWarning)
    line.collimators.assign_optics(nemitt_x=nemitt_x, nemitt_y=nemitt_y, twiss=twiss)

def get_optics_at(names, *, twiss=None, line=None):
    warn("The function xcoll.get_optics_at() is deprecated and will be "
       + "removed in the future. Please use line.collimators.get_optics_at() instead.", FutureWarning)
    return line.collimators.get_optics_at(names=names, twiss=twiss)

def open_collimators(line, names=None):
    warn("The function xcoll.open_collimators() is deprecated and will be "
       + "removed in the future. Please use line.collimators.open_collimators() instead.", FutureWarning)
    line.collimators.open(names=names)

def send_to_parking(line, names=None):
    warn("The function xcoll.send_to_parking() is deprecated and will be "
       + "removed in the future. Please use line.collimators.send_to_parking() instead.", FutureWarning)
    line.collimators.to_parking(names=names)

def enable_scattering(line):
    warn("The function xcoll.enable_scattering() is deprecated and will be "
       + "removed in the future. Please use line.scattering.enable() instead.", FutureWarning)
    line.scattering.enable()

def disable_scattering(line):
    warn("The function xcoll.disable_scattering() is deprecated and will be "
       + "removed in the future. Please use line.scattering.disable() instead.", FutureWarning)
    line.scattering.disable()
