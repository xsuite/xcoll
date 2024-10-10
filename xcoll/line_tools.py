# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import warnings
warnings.simplefilter("always")

import xtrack as xt

from .beam_elements import element_classes, _all_collimator_classes


class XcollCollimatorAPI:
    def __init__(self, line):
        self._line = line

    @property
    def line(self):
        return self._line


    def get_optics_at(self, names, *, twiss=None):
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

    def assign_optics(self, nemitt_x=None, nemitt_y=None, twiss=None):
        if not self.line._has_valid_tracker():
            raise Exception("Please build tracker before setting the openings!")
        names = self.line.get_elements_of_type(_all_collimator_classes)[1]
        tw_upstream, tw_downstream = self.get_optics_at(names, twiss=twiss)
        beta_gamma_rel = self.line.particle_ref._xobject.gamma0[0]*self.line.particle_ref._xobject.beta0[0]
        for coll in names:
            self.line[coll].assign_optics(name=coll, nemitt_x=nemitt_x, nemitt_y=nemitt_x, twiss_upstream=tw_upstream,
                                    twiss_downstream=tw_downstream, beta_gamma_rel=beta_gamma_rel)

    def open(self, names=None):
        if names is None:
            names = self.line.get_elements_of_type(_all_collimator_classes)[1]
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




def assign_optics_to_collimators(line, nemitt_x=None, nemitt_y=None, twiss=None):
    warnings.warn("The function xcoll.assign_optics_to_collimators() is deprecated and will be "
                + "removed in the future. Please use line.scattering.assign_optics() instead.", DeprecationWarning)
    line.collimators.assign_optics(nemitt_x=nemitt_x, nemitt_y=nemitt_y, twiss=twiss)

def get_optics_at(names, *, twiss=None, line=None):
    warnings.warn("The function xcoll.get_optics_at() is deprecated and will be "
                + "removed in the future. Please use line.scattering.get_optics_at() instead.", DeprecationWarning)
    return line.collimators.get_optics_at(names=names, twiss=twiss)

def open_collimators(line, names=None):
    warnings.warn("The function xcoll.open_collimators() is deprecated and will be "
                + "removed in the future. Please use line.scattering.open_collimators() instead.", DeprecationWarning)
    line.collimators.open(names=names)

def send_to_parking(line, names=None):
    warnings.warn("The function xcoll.send_to_parking() is deprecated and will be "
                + "removed in the future. Please use line.scattering.send_to_parking() instead.", DeprecationWarning)
    line.collimators.to_parking(names=names)

def enable_scattering(line):
    warnings.warn("The function xcoll.enable_scattering() is deprecated and will be "
                + "removed in the future. Please use line.scattering.enable() instead.", DeprecationWarning)
    line.scattering.enable()

def disable_scattering(line):
    warnings.warn("The function xcoll.disable_scattering() is deprecated and will be "
                + "removed in the future. Please use line.scattering.disable() instead.", DeprecationWarning)
    line.scattering.disable()
