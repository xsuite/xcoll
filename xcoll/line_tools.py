# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import xtrack as xt

from .beam_elements import element_classes, _all_collimator_classes


def assign_optics_to_collimators(line, nemitt_x=None, nemitt_y=None, twiss=None):
    if not line._has_valid_tracker():
        raise Exception("Please build tracker before setting the openings!")
    names = line.get_elements_of_type(_all_collimator_classes)[1]
    tw_upstream, tw_downstream = get_optics_at(names, twiss=twiss, line=line)
    beta_gamma_rel = line.particle_ref._xobject.gamma0[0]*line.particle_ref._xobject.beta0[0]
    for coll in names:
        line[coll].assign_optics(name=coll, nemitt_x=nemitt_x, nemitt_y=nemitt_x, twiss_upstream=tw_upstream,
                                 twiss_downstream=tw_downstream, beta_gamma_rel=beta_gamma_rel)

def get_optics_at(names, *, twiss=None, line=None):
    if twiss is None:
        if not line._has_valid_tracker():
            raise Exception("Please pass a line and build tracker before computing the optics for the openings!")
        twiss = line.twiss()
    if not hasattr(names, '__iter__') and not isinstance(names, str):
        names = [names]
    coll_entry_mask = twiss.mask[names]
    tw_entry = twiss.rows[coll_entry_mask]
    tw_exit = twiss.rows[coll_entry_mask+1]
    tw_exit.name = tw_entry.name
    return tw_entry, tw_exit


def open_collimators(line, names=None):
    if names is None:
        names = line.get_elements_of_type(_all_collimator_classes)[1]
    if len(names) == 0:
        print("No collimators found in line.")
    else:
        for coll in names:
            line[coll].open_jaws(keep_tilts=False)
            line[coll].gap = None

def send_to_parking(line, names=None):
    if names is None:
        names = line.get_elements_of_type(_all_collimator_classes)[1]
    if len(names) == 0:
        print("No collimators found in line.")
    else:
        raise NotImplementedError("Need to move this to new type manager or so.")


def enable_scattering(line):
    elements = line.get_elements_of_type(element_classes)[0]
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
            el.enable_scattering()
        # self.line.tracker.io_buffer = self._io_buffer
        # self._set_record_interaction_record()

def disable_scattering(line):
    elements = line.get_elements_of_type(element_classes)[0]
    if len(elements) == 0:
        print("No xcoll elements found in line.")
    else:
        for el in elements:
            el.disable_scattering()
