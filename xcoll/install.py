# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from warnings import warn

# This file is deprecated and will be removed in the future.


def install_elements(line, names, elements, *, at_s=None, apertures=None, need_apertures=False, s_tol=1.e-6):
    warn("The function xcoll.install_elements() is deprecated and will be removed in the future. "
       + "Please use line.collimators.install() instead.", FutureWarning)
    line.collimators.install(names=names, elements=elements, at_s=at_s, apertures=apertures,
                             need_apertures=need_apertures, s_tol=s_tol)

def check_element_position(line, name, s_start, s_end, at_s, length=None, s_tol=1.e-6):
    warn("The function xcoll.check_element_position() is deprecated and will be removed "
       + "in the future. Please use line.collimators.check_position() instead.", FutureWarning)
    return line.collimators.check_position(name, s_start, s_end, at_s, length=length, s_tol=s_tol)

def get_apertures_at_s(tab, s, s_tol=1.e-6):
    warn("The function xcoll.get_apertures_at_s() is deprecated and will be removed in the future. "
       + "Please use line.collimators.get_apertures_at_s() instead.", FutureWarning)
    return get_apertures_at_s(s, s_tol=s_tol)

def get_aperture_for_element(line, name, s_start, s_end, aperture=None, tab=None, s_tol=1.e-6):
    warn("The function xcoll.get_aperture_for_element() is deprecated and will be removed in the future. "
       + "Please use line.collimators.get_aperture() instead.", FutureWarning)
    return line.get_aperture(name, s_start, s_end, aperture=aperture, tab=tab, s_tol=s_tol)

def prepare_space_for_element(line, name, s_start, s_end, tab=None, s_tol=1.e-6):
    warn("The function xcoll.prepare_space_for_element() is deprecated and will be removed in the future. "
       + "Please use line.collimators.prepare_space() instead.", FutureWarning)
    return line.collimators.prepare_space(name, s_start, s_end, tab=tab, s_tol=s_tol)
