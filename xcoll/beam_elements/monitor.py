# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from ..general import _pkg_root

class EmittanceMonitorRecord(xo.Struct):
    count      = xo.Float64[:]
    x_sum1     = xo.Float64[:]
    px_sum1    = xo.Float64[:]
    y_sum1     = xo.Float64[:]
    py_sum1    = xo.Float64[:]
    x_x_sum2   = xo.Float64[:]
    x_px_sum2  = xo.Float64[:]
    x_y_sum2   = xo.Float64[:]
    x_py_sum2  = xo.Float64[:]
    px_px_sum2 = xo.Float64[:]
    px_y_sum2  = xo.Float64[:]
    px_py_sum2 = xo.Float64[:]
    y_y_sum2   = xo.Float64[:]
    y_py_sum2  = xo.Float64[:]
    py_py_sum2 = xo.Float64[:]

class EmittanceMonitor(xt.BeamElement):
    _xofields={
        'part_id_start':      xo.Int64,
        'part_id_end':        xo.Int64,
        'start_at_turn':      xo.Int64,
        'stop_at_turn':       xo.Int64,
        'frev':               xo.Float64,
        'sampling_frequency': xo.Float64,
        '_index':             xt.RecordIndex,
        'data':               EmittanceMonitorRecord,
        '_cached':             xo.Int8
    }

    behaves_like_drift = True
    allow_loss_refinement = True

    _extra_c_sources = [
        xt._pkg_root.joinpath('headers/atomicadd.h'),
        _pkg_root.joinpath('beam_elements/elements_src/emittance_monitor.h')
    ]

    def __init__(self, **kwargs):
        """
        Monitor to save the normalised beam emittance

        Similar to the BeamSizeMonitor and BeamPositionMonitor, it allows for
        arbitrary sampling rate and can thus not only be used to monitor bunch
        emittance, but also to record coasting beams. See their documentation
        for more information on how to use `frev` and `sampling_frequency`.

        Args:
            num_particles (int, optional): Number of particles to monitor,
                starting from 0. Defaults to -1 which means ALL.
            particle_id_range (tuple, optional): Range of particle ids to
                monitor (start, stop). Stop is exclusive. Defaults to
                (particle_id_start, particle_id_start+num_particles).
            start_at_turn (int): First turn of reference particle (inclusive)
                at which to monitor. Defaults to 0.
            stop_at_turn (int): Last turn of reference particle (exclusive) at
                which to monitor. Defaults to start_at_turn + 1.
            frev (float): Revolution frequency in Hz of circulating beam (used
                to relate turn number to sample index). Defaults to 1.
            sampling_frequency (float): Sampling frequency in Hz. Defaults to 1.
        """
        if '_xobject' not in kwargs:
            if 'particle_id_range' in kwargs:
                assert 'num_particles' not in kwargs
                particle_id_range = kwargs.pop('particle_id_range')
                kwargs['part_id_start'] = int(particle_id_range[0])
                kwargs['part_id_end'] = int(particle_id_range[1])
            else:
                num_particles = int(kwargs.pop('num_particles', -1))
                if num_particles == -1:
                    kwargs['part_id_start'] = 0
                    kwargs['part_id_end'] = -1
                else:
                    kwargs['part_id_start'] = kwargs.pop('particle_id_start', 0)
                    kwargs['part_id_end'] = kwargs['part_id_start'] + num_particles
            kwargs['start_at_turn'] = int(kwargs.get('start_at_turn', 0))
            kwargs['stop_at_turn']  = int(kwargs.get('stop_at_turn', kwargs['start_at_turn']+1))
            kwargs.setdefault('frev', 1.)
            kwargs.setdefault('sampling_frequency', 1.)
            kwargs['_cached'] = False
            if "data" not in kwargs:
                # explicitely init with zeros (instead of size only) to have consistent initial values
                size = int(round((kwargs['stop_at_turn'] - kwargs['start_at_turn']) \
                                  * kwargs['sampling_frequency'] / kwargs['frev']))
                kwargs['data'] = {field.name: np.zeros(size) for field in EmittanceMonitorRecord._fields}
        super().__init__(**kwargs)


    @property
    def emitt_x(self):
        if not self._cached:
            self._calculate_all()
        return self._emitt_x

    @property
    def emitt_y(self):
        if not self._cached:
            self._calculate_all()
        return self._emitt_y

    @property
    def nemitt_x(self):
        if not self._cached:
            self._calculate_all()
        if not hasattr(self, '_beta_gamma_rel'):
            raise ValueError("Need to call `set_beta_gamma_rel()` first!")
        return self._emitt_x * self._beta_gamma_rel

    @property
    def nemitt_y(self):
        if not self._cached:
            self._calculate_all()
        if not hasattr(self, '_beta_gamma_rel'):
            raise ValueError("Need to call `set_beta_gamma_rel()` first!")
        return self._emitt_y * self._beta_gamma_rel


    def set_beta_gamma_rel(self, particles=None, beta=None, gamma=None):
        if particles is not None:
            if beta is not None or gamma is not None:
                raise ValueError("Use either `particles` or `beta` and `gamma`!")
            beta = particles.beta0[0]
            gamma = particles.gamma0[0]
        self._beta_gamma_rel = beta * gamma


    def _calculate_all(self):
        self._cached = True

        # Calculate mean, variance, and std
        N = self.count
        with np.errstate(invalid='ignore'):  # NaN for zero particles is expected behaviour
            for field in [f.name for f in EmittanceMonitorRecord._fields]:
                if field.endswith('_sum1'):
                    x = field[:-5]
                    mean = getattr(self, field) / N
                    setattr(self, f'_{x}_mean', mean)
                elif field.endswith('_sum2'):
                    x1, x2 = field[:-5].split('_')
                    mean1 = getattr(self, f'{x1}_sum1') / N
                    mean2 = getattr(self, f'{x2}_sum1') / N
                    sum2 = getattr(self, field)
                    variance = sum2 / (N - 1) - mean1 * mean2 * N / (N - 1)
                    setattr(self, f'_{x1}_{x2}_var', variance)
                    if x1 == x2:
                        setattr(self, f'_{x1}_std', np.sqrt(variance))

        # Calculate emittances
        emitt_x = np.sqrt(self.x_x_var * self.px_px_var - self.x_px_var**2)
        emitt_y = np.sqrt(self.y_y_var * self.py_py_var - self.y_py_var**2)
        setattr(self, '_emitt_x', emitt_x)
        setattr(self, '_emitt_y', emitt_y)


    def __getattr__(self, attr):
        if attr in [f.name for f in EmittanceMonitorRecord._fields]:
            return getattr(self.data, attr).to_nparray()

        else:
            if attr.startswith('_'):
                raise AttributeError(f"Attribute {attr} not set!")

            if attr.endswith('_mean') or attr.endswith('_var') or attr.endswith('_std'):
                if not self._cached:
                    self._calculate_all()
                return getattr(self, f'_{attr}')
            else:
                raise AttributeError(f"EmittanceMonitor has no attribute '{attr}'")
