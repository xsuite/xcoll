# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import xtrack as xt

from .beam_elements import element_classes

def install_elements(line, names, elements, *, at_s=None, apertures=None, need_apertures=False, s_tol=1.e-6):
    if line._has_valid_tracker():
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
        assert isinstance(el, element_classes)
        el._tracking = False

    # Get positions
    tab = line.get_table()
    tt = tab.rows[[name for name in names if name in line.element_names]]
    s_start = []
    for name, s, l in zip(names, at_s, length):
        if s is None:
            s_start.append(_get_s_start(line, name, l, tt))
        else:
            s_start.append(s)
    s_start = np.array(s_start)
    s_end = s_start + length

    # Check positions
    l_line = line.get_length() 
    for s1, s2, name, s3 in zip(s_start, s_end, names, at_s):
        check_element_position(line, name, s1, s2, s3, l_line, s_tol=s_tol)

    # Look for apertures
    aper_upstream   = []
    aper_downstream = []
    for s1, s2, name, aper in zip(s_start, s_end, names, apertures):
        if not need_apertures:
            aper_upstream.append(None)
            aper_downstream.append(None)
        else:
            aper1, aper2 = get_aperture_for_element(line, name, s1, s2, aper, tab, s_tol=s_tol)
            aper_upstream.append(aper1)
            aper_downstream.append(aper2)

    # Remove elements at location of collimator (by changing them into markers)
    for s1, s2, name in zip(s_start, s_end, names):
        prepare_space_for_element(line, name, s1, s2, tab=tab, s_tol=s_tol)

    # Install
    line._insert_thick_elements_at_s(element_names=list(names), elements=elements, at_s=s_start, s_tol=s_tol)

    # Install apertures
    if need_apertures:
        for s1, name, aper1, aper2 in zip(s_start, names, aper_upstream, aper_downstream):
            line.insert_element(element=aper1, name=f'{name}_aper_upstream', at=name, s_tol=s_tol)
            idx = line.element_names.index(name) + 1
            line.insert_element(element=aper2, name=f'{name}_aper_downstream', at=idx, s_tol=s_tol)


def _get_s_start(line, name, length, tab=None):
    if tab is None:
        tab = line.get_table()
    if name in line.element_names and hasattr(line[name], 'length'):
        existing_length = line[name].length
    else:
        existing_length = 0
    return tab.rows[name].s[0] + existing_length/2. - length/2


def check_element_position(line, name, s_start, s_end, at_s, length=None, s_tol=1.e-6):
    if at_s is None:
        if name not in line.element_names:
            raise ValueError(f"Element {name} not found in line. Provide `at_s`.")
    elif name in line.element_names:
        if at_s < s_start or at_s > s_end:
            raise ValueError(f"Element {name} already exists in line at different "
                           + f"location: at_s = {at_s}, exists at [{s_start}, {s_end}].")
    if length is None:
        length = line.get_length() 
    if s_start <= s_tol:
        raise ValueError(f"Position of {name} too close to start of line. Please cycle.")
    if s_end >= length - s_tol:
        raise ValueError(f"Position of {name} too close to end of line. Please cycle.")


def get_apertures_at_s(tab, s, s_tol=1.e-6):
    tab_s = tab.rows[s-s_tol:s+s_tol:'s']
    aper = tab_s.rows[[cls.startswith('Limit') for cls in tab_s.element_type]]
    if len(aper) == 0:
        return None
    elif len(aper) == 1:
        return aper.name[0]
    else:
        raise ValueError(f"Multiple apertures found at location {s} with "
                       + f"tolerance {s_tol}: {aper.name}. Not supported.")


def get_aperture_for_element(line, name, s_start, s_end, aperture=None, tab=None, s_tol=1.e-6):
    if aperture is not None:
        if isinstance(aperture, str):
            aper1 = line[aperture]
            aper2 = line[aperture]
        elif hasattr(aperture, '__iter__'):
            if len(aperture) != 2:
                raise ValueError(f"The value `aperture` should be None or a list "
                               + f"[upstream, downstream].")
            assert aperture[0] is not None and aperture[1] is not None
            if isinstance(aperture[0], str):
                aper1 = line[aperture[0]]
            if isinstance(aperture[1], str):
                aper2 = line[aperture[1]]
        else:
            aper1 = aperture
            aper2 = aperture
        if not xt.line._is_aperture(aper1, line):
            raise ValueError(f"Not a valid aperture: {aper1}")
        if not xt.line._is_aperture(aper2, line):
            raise ValueError(f"Not a valid aperture: {aper2}")
        return aper1.copy(), aper2.copy()
    else:
        if tab is None:
            tab = line.get_table()
        aper1 = get_apertures_at_s(tab, s_start, s_tol=s_tol)
        aper2 = get_apertures_at_s(tab, s_end, s_tol=s_tol)
        if aper1 is None and aper2 is not None:
            aper1 = aper2
            print(f"Warning: Could not find upstream aperture for {name}! "
                + f"Used copy of downstream aperture. Proceed with caution.")
        elif aper2 is None and aper1 is not None:
            aper2 = aper1
            print(f"Warning: Could not find downstream aperture for {name}! "
                + f"Used copy of upstream aperture. Proceed with caution.")
        elif aper1 is None and aper2 is None:
            aper_mid = get_apertures_at_s(tab, (s_start+s_end)/2, s_tol=s_tol)
            if aper_mid is None:
                raise ValueError(f"No aperture found for {name}! Please provide one.")
            if line[aper_mid].allow_rot_and_shift \
            and xt.base_element._tranformations_active(line[aper_mid]):
                print(f"Warning: Using the centre aperture for {name}, but "
                    + f"transformations are present. Proceed with caution.")
            aper1 = aper_mid
            aper2 = aper_mid
        return line[aper1].copy(), line[aper2].copy()


def prepare_space_for_element(line, name, s_start, s_end, tab=None, s_tol=1.e-6):
    if tab is None:
        tab = line.get_table()
    tt = tab.rows[s_start-s_tol:s_end+s_tol:'s']
    for element_name, element_type in zip(tt.name[:-1], tt.element_type[:-1]):
        if element_type == 'Marker' or element_type.startswith('Drift'):
            continue
        if not element_type.startswith('Limit'):
            print(f"Warning: Removed active element {element_name} "
                + f"at location inside collimator!")
        length = line[element_name].length if hasattr(line[element_name], 'length') else 0
        line.element_dict[element_name] = xt.Drift(length=length)

