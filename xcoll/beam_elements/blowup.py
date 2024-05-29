# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..general import _pkg_root
from .base import InvalidXcoll


class BlowUp(InvalidXcoll):
    _xofields = {
        '_plane':       xo.Int8,
        '_amplitude':   xo.Float64,  # This is the one used in C, it is already de-calibrated
        '_calibration': xo.Float64,
        '_active':      xo.Int8,
    }

    isthick = False
    behaves_like_drift = False
    allow_track = True
    needs_rng = True
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _skip_in_to_dict  = ['_amplitude', '_plane', '_calibration', '_active']
    _store_in_to_dict = ['amplitude', 'plane', 'calibration']

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','blowup.h')
    ]

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            if 'plane' not in kwargs:
                raise ValueError("The plane of the BlowUp must be specified.")
            to_assign['plane']       = kwargs.pop('plane')
            to_assign['calibration'] = kwargs.pop('calibration', 1.)
            kwargs['_calibration'] = 1.
            if 'amplitude' not in kwargs:
                raise ValueError("The amplitude of the BlowUp must be specified "
                               + "(in terms of calibrated units.")
            to_assign['amplitude']   = kwargs.pop('amplitude')
            kwargs['_active'] = 0
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)
        self._name = None

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return self.__class__(length=-self.length,
                              _context=_context, _buffer=_buffer, _offset=_offset)

    @property
    def plane(self):
        if self._plane == 1:
            return 'H'
        elif self._plane == -1:
            return 'V'
        else:
            raise ValueError("The plane of the BlowUp is not recognized.")

    @plane.setter
    def plane(self, val):
        if val.lower() == 'h':
            self._plane = 1
        elif val.lower() == 'v':
            self._plane = -1
        else:
            raise ValueError("The plane of the BlowUp must be either 'H' or 'V'.")

    @property
    def amplitude(self):
        return self._amplitude / self._calibration

    @amplitude.setter
    def amplitude(self, val):
        if val< 0:
            raise ValueError("The amplitude cannot be negative!")
        self._amplitude = val * self._calibration

    @property
    def calibration(self):
        return self._calibration

    @calibration.setter
    def calibration(self, val):
        if val< 0:
            raise ValueError("The calibration cannot be negative!")
        previous_amplitude = self.amplitude
        self._calibration = val
        self.amplitude = previous_amplitude

    def calibrate(self, twiss, nemitt_x, nemitt_y, name=None):
        if name is None:
            if self._name is None:
                raise ValueError("The name of the BlowUp must be provided.")
            name = self._name
        beam_sizes = twiss.get_beam_covariance(nemitt_x=nemitt_x, nemitt_y=nemitt_y)
        sigma_p = beam_sizes.rows[name][f"sigma_p{'x' if self.plane == 'H' else 'y'}"][0]
        self.calibration = sigma_p*5./np.sqrt(500)  # Gain 4-5 sigma average in 500 turns

    def activate(self):
        self._active = True

    def deactivate(self):
        self._active = False

    def install(self, line, name, *, at_s=None, at=None, need_apertures=True, aperture=None, s_tol=1.e-6):
        line.insert_element(element=self, name=name, at_s=at_s, at=at, s_tol=s_tol)
        self._name = name
        if need_apertures:
            if aperture is not None:
                aper_upstream   = aperture.copy()
                aper_downstream = aperture.copy()
            else:
                idx = line.element_names.index(name)
                while True:
                    if xt.line._is_aperture(line.elements[idx], line):
                        aper_upstream = line.elements[idx].copy()
                        break
                    idx -= 1
                idx = line.element_names.index(name)
                while True:
                    if xt.line._is_aperture(line.elements[idx], line):
                        aper_downstream = line.elements[idx].copy()
                        break
                    idx += 1
            line.insert_element(element=aper_upstream, name=f'{name}_aper_upstream', at=name, s_tol=s_tol)
            idx = line.element_names.index(name) + 1
            line.insert_element(element=aper_downstream, name=f'{name}_aper_downstream', at=idx, s_tol=s_tol)
