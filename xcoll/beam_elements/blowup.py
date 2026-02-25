# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt


class BlowUp(xt.BeamElement):
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

    _isthick = False
    needs_rng = True
    prototype = "blowup"
    allow_track = True
    iscollective = False
    has_backtrack = False
    behaves_like_drift = False
    allow_rot_and_shift = True
    allow_loss_refinement = False
    name_associated_aperture = None
    skip_in_loss_location_refinement = False

    _noexpr_fields = {'plane', 'name', 'line'}
    _skip_in_to_dict = ['_max_kick', '_plane', '_calibration', '_active']
    _store_in_to_dict = ['amplitude', 'plane', 'calibration', 'beta0',
                         'gamma0', 'name']

    _depends_on = [xt.RandomUniform]
    _extra_c_sources = [
        '#include "xcoll/beam_elements/elements_src/blowup.h"'
    ]


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            start_at_turn = int(kwargs.get('start_at_turn', 0))
            stop_at_turn  = int(kwargs.get('stop_at_turn', start_at_turn+1))
            kwargs['start_at_turn'] = start_at_turn
            kwargs['stop_at_turn'] = stop_at_turn
            rans = np.random.uniform(size=stop_at_turn-start_at_turn)
            kwargs['_rans'] = 2*rans - 1
            if 'plane' in kwargs:
                to_assign['plane'] = kwargs.pop('plane')
            to_assign['calibration'] = kwargs.pop('calibration', 1.)
            to_assign['amplitude'] = kwargs.pop('amplitude', 1)
            kwargs['_calibration'] = 1.
            to_assign['_beta0'] = kwargs.pop('beta0', None)
            to_assign['_gamma0'] = kwargs.pop('gamma0', None)
            to_assign['_name'] = kwargs.pop('name', None)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


    @classmethod
    def install(cls, line, name, *, at_s=None, at=None, need_apertures=True,
                aperture=None, s_tol=1.e-6, **kwargs):
        """Shortcut to install a BlowUp in a line, which also sets the
        optics parameters. If `need_apertures` is True, the BlowUp will
        be installed together with upstream and downstream apertures.
        """
        self = cls(**kwargs)
        if name in line.element_names:
            raise ValueError(f"Element {name} already exists in the line as "
                             f"{line[name].__class__.__name__}.")
        line.insert_element(element=self, name=name, at_s=at_s, at=at,
                            s_tol=s_tol)
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
            line.insert_element(element=aper_upstream, at=name, s_tol=s_tol,
                                name=f'{name}_aper_upstream')
            idx = line.element_names.index(name) + 1
            line.insert_element(element=aper_downstream, at=idx, s_tol=s_tol,
                                name=f'{name}_aper_downstream')
        if hasattr(line, 'particle_ref') and line.particle_ref is not None:
            self._beta0 = line.particle_ref.beta0[0]
            self._gamma0 = line.particle_ref.gamma0[0]
        return self


    @property
    def beta0(self):
        return self._beta0

    @property
    def gamma0(self):
        return self._gamma0

    @property
    def name(self):
        return self._name


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
            raise ValueError("The plane of the BlowUp must be either 'H' or "
                             "'V'.")

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


    def calibrate_by_emittance(self, nemitt, *, line=None, twiss=None,
                               name=None, beta0=None, gamma0=None):
        # The emittance gained per turn of blow up can be calculated as
        # follows:
        #     <px^2> = gamma * eps
        #     gamma * eps = <px^2> = <(px0 + dpx)^2>
        # where px0 is the initial momentum and dpx is the kick, and
        # eps0 and eps are the geometric emittance before resp. after
        # the kick, and gamma is the Courant-Snyder parameter. If each
        # particle gets another random kick, we can write:
        #     <(px0 + dpx)^2> = <px0^2> + 2<px0 dpx> + <dpx^2>
        #                     = gamma * eps0 + 2<px0 dpx> + 1/3 dpx_max^2    (variance of uniform distribution is 1/12(b-a)^2)
        # The covariance can be estimated with an upper bound:
        #     0  <=  2<px0 dpx>  <=  2*sqrt(<px0^2> <dpx^2>)
        #                         = 2/sqrt(3) * sqrt(gamma * eps0) * dpx_max
        # So that we finally get for the emittance growth after one turn
        # of blow-up:
        #     deps = (eps-eps0)/eps0 = [dPM^2, dPM^2 + 2dPM]    (dPM = dpx_max/sqrt(3 * gamma * eps0))
        # Or, to get an emittance growth of deps, we need:
        #     dpx_max = [sqrt(3* gamma * eps0 * deps), sqrt(3 * gamma * eps0)*(sqrt(1+deps) -1)]
        # If all particles get the same kick, we can write:
        #     <(px0 + dpx)^2> = <px0^2> + 2<px0> dpx
        # However, now things get a bit more complicated. At the start
        # <px0> is the closed orbit, but when the bunch has had a kick
        # in the previous turn, the bunch is no longer centred around
        # the closed orbit. We can write:
        #     deps = (eps-eps0)/eps0 = 2<px0> dpx/gamma/eps0
        if beta0 is None or gamma0 is None:
            if self._beta0 is not None and self._gamma0 is not None:
                beta0 = self._beta0
                gamma0 = self._gamma0
            elif line is not None and hasattr(line, 'particle_ref') \
            and line.particle_ref is not None:
                beta0 = line.particle_ref.beta0[0]
                gamma0 = line.particle_ref.gamma0[0]
            else:
                raise ValueError("Either a line with a particle_ref, or "
                                 "beta0 and gamma0 must be provided!")
        if twiss is None:
            if line is None:
                raise ValueError("Either `line` or `twiss` must be provided!")
            twiss = line.twiss(reverse=False)
        coord = f"gam{'x' if self.plane == 'H' else 'y'}"
        if name is None:
            name = self.name
        gamma = twiss.rows[name][coord][0]
        # Asumming a Gaussian beam, we will have shifted all beam beyond
        # 5 sigma with an emittance growth of around a factor 60. In
        # practice, a ~50x increase is enough due to the covariance
        # terms. We want to do this slowly, over 1000 turns.
        deps = 50/1000
        gemitt = nemitt / beta0 / gamma0
        self.calibration = np.sqrt(3*gamma*gemitt*deps)

    def activate(self):
        self._active = True

    def deactivate(self):
        self._active = False
