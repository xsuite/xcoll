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
        '_kick_rms':    xo.Float64,
        '_calibration': xo.Float64,
        '_active':      xo.Int8,
    }

    isthick = False
    behaves_like_drift = False
    allow_track = True
    needs_rng = True
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _skip_in_to_dict  = ['_kick_rms', '_plane', '_calibration', '_active']
    _store_in_to_dict = ['amplitude', 'plane', 'calibration']

    _depends_on = [InvalidXcoll, xt.RandomNormal]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','blowup.h')
    ]

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            if 'plane' in kwargs:
                to_assign['plane']   = kwargs.pop('plane')
            to_assign['calibration'] = kwargs.pop('calibration', 1.)
            to_assign['amplitude']   = kwargs.pop('amplitude', 1)
            kwargs['_calibration']   = 1.
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)
        if not hasattr(self, '_name'):
            self._name = None
        if not hasattr(self, '_line'):
            self._line = None

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
        return (self._kick_rms / self.calibration)**2

    @amplitude.setter
    def amplitude(self, val):
        if val< 0:
            raise ValueError("The amplitude cannot be negative!")
        self._kick_rms = np.sqrt(val) * self.calibration

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

    def calibrate_by_emittance(self, nemitt, *, name=None, line=None,
                               twiss=None, beta_gamma_rel=None):
        # Emittance grows as < emittance * dpx / px>
        # We approximate this as <emittance> * <dpx> / sqrt(<px^2>), and using
        # <px^2> = gamma * emittance_geom
        if name is None:
            if self._name is None:
                raise ValueError("The name of the BlowUp must be provided.")
            name = self._name
        if line is None:
            line = self._line
        if beta_gamma_rel is None:
            if line is None:
                raise ValueError("Need to provide `beta_gamma_el` or a line with a reference particle.")
            elif not hasattr(line, 'particle_ref'):
                raise ValueError("The provided line has no reference particle.")
            beta_gamma_rel = line.particle_ref._xobject.gamma0[0]*line.particle_ref._xobject.beta0[0]
        if twiss is None:
            if line is None:
                raise ValueError("Need to provide `twiss` or a line.")
            twiss = line.twiss()
        gamma = twiss.rows[name][f"gam{'x' if self.plane == 'H' else 'y'}"][0]
        # Small correction of 20% (empirical) to account for underestimation due to statistics error
        self.calibration = np.sqrt(gamma*nemitt/beta_gamma_rel) / 1.2

    def activate(self):
        self._active = True

    def deactivate(self):
        self._active = False

    def install(self, line, name, *, at_s=None, at=None, need_apertures=True, aperture=None, s_tol=1.e-6):
        line.insert_element(element=self, name=name, at_s=at_s, at=at, s_tol=s_tol)
        self._name = name
        self._line = line
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
