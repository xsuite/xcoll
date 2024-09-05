# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..general import _pkg_root
from .base import InvalidXcoll


class BlowUp(InvalidXcoll):
    _xofields = {
        '_plane':               xo.Int8,
        'start_at_turn':        xo.Int64,
        'stop_at_turn':         xo.Int64,
        'use_individual_kicks': xo.Int8,
        '_max_kick':            xo.Float64,
        '_rans':                xo.Float64[:],
        '_calibration':         xo.Float64,
        '_active':              xo.Int8,
    }

    isthick = False
    behaves_like_drift = False
    allow_track = True
    needs_rng = True
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _skip_in_to_dict  = ['_max_kick', '_plane', '_calibration', '_active']
    _store_in_to_dict = ['amplitude', 'plane', 'calibration']

    _depends_on = [InvalidXcoll, xt.RandomUniform]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','blowup.h')
    ]

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            start_at_turn = int(kwargs.get('start_at_turn', 0))
            stop_at_turn  = int(kwargs.get('stop_at_turn', start_at_turn+1))
            kwargs['start_at_turn'] = start_at_turn
            kwargs['stop_at_turn']  = stop_at_turn
            kwargs['_rans'] = 2*np.random.uniform(size=stop_at_turn-start_at_turn) - 1
            if 'plane' in kwargs:
                to_assign['plane']   = kwargs.pop('plane')
            to_assign['calibration'] = kwargs.pop('calibration', 1.)
            to_assign['amplitude']   = kwargs.pop('amplitude', 1)
            kwargs['_calibration']   = 1.
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


    @classmethod
    def install(cls, line, name, *, at_s=None, at=None, need_apertures=True, aperture=None, s_tol=1.e-6, **kwargs):
        self = cls(**kwargs)
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
        return self

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return xt.Marker(_context=_context, _buffer=_buffer, _offset=_offset)


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
        return (self._max_kick / self.calibration)**2

    @amplitude.setter
    def amplitude(self, val):
        if val< 0:
            raise ValueError("The amplitude cannot be negative!")
        self._max_kick = np.sqrt(val) * self.calibration

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

    @property
    def name(self):
        if not hasattr(self, '_name'):
            raise ValueError("Name not set! Install the blow-up using xc.BlowUp.install() "
                             "or manually set the name after installation.")
        return self._name

    @name.setter
    def name(self, val):
        self._name = val

    @property
    def line(self):
        if not hasattr(self, '_line'):
            raise ValueError("Line not set! Install the blow-up using xc.BlowUp.install() "
                             "or manually set the line after installation.")
        return self._line

    @line.setter
    def line(self, val):
        self._line = val


    def calibrate_by_emittance(self, nemitt, *, twiss=None, beta_gamma_rel=None):
        # The emittance gained per turn of blow up can be calculated as follows:
        #     <px^2> = gamma * eps
        #     gamma * eps = <px^2> = <(px0 + dpx)^2>
        # where px0 is the initial momentum and dpx is the kick, and eps0 and eps
        # are the geometric emittance before resp. after the kick, and gamma is the
        # Courant-Snyder parameter.
        # If each particle gets another random kick, we can write:
        #     <(px0 + dpx)^2> = <px0^2> + 2<px0 dpx> + <dpx^2>
        #                     = gamma * eps0 + 2<px0 dpx> + 1/3 dpx_max^2    (variance of uniform distribution is 1/12(b-a)^2)
        # The covariance can be estimated with an upper bound:
        #     0  <=  2<px0 dpx>  <=  2*sqrt(<px0^2> <dpx^2>)
        #                         = 2/sqrt(3) * sqrt(gamma * eps0) * dpx_max
        # So that we finally get for the emittance growth after one turn of blow-up:
        #     deps = (eps-eps0)/eps0 = [dPM^2, dPM^2 + 2dPM]    (dPM = dpx_max/sqrt(3 * gamma * eps0))
        # Or, to get an emittance growth of deps, we need:
        #     dpx_max = [sqrt(3* gamma * eps0 * deps), sqrt(3 * gamma * eps0)*(sqrt(1+deps) -1)]
        # If all particles get the same kick, we can write:
        #     <(px0 + dpx)^2> = <px0^2> + 2<px0> dpx
        # However, now things get a bit more complicated. At the start <px0> is the closed orbit,
        # but when the bunch has had a kick in the previous turn, the bunch is no longer centred
        # around the closed orbit. We can write:
        #     deps = (eps-eps0)/eps0 = 2<px0> dpx/gamma/eps0
        name = self.name
        line = self.line
        if beta_gamma_rel is None:
            if not hasattr(line, 'particle_ref'):
                raise ValueError("The provided line has no reference particle. Use the argument `beta_gamma_rel`")
            beta_gamma_rel = line.particle_ref.gamma0[0]*line.particle_ref.beta0[0]
        if twiss is None:
            twiss = line.twiss()
        gamma = twiss.rows[name][f"gam{'x' if self.plane == 'H' else 'y'}"][0]
        # Asumming a Gaussian beam, we will have shifted all beam beyond 5 sigma with an emittance growth
        # of around a factor 60. In practive, a ~50x increase is enough due to the covariance terms.
        # We want to do this slowly, over 1000 turns.
        deps = 50/1000
        self.calibration = np.sqrt(3*gamma*nemitt/beta_gamma_rel*deps)

    def activate(self):
        self._active = True

    def deactivate(self):
        self._active = False
