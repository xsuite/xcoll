# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import io
import numpy as np

import xobjects as xo
import xtrack as xt

from .. import json
from ..compare import deep_equal


class ParticleStatsMonitorRecord(xo.Struct):
    count            = xo.Int64[:]
    x_sum1           = xo.Float64[:]
    px_sum1          = xo.Float64[:]
    y_sum1           = xo.Float64[:]
    py_sum1          = xo.Float64[:]
    zeta_sum1        = xo.Float64[:]
    pzeta_sum1       = xo.Float64[:]
    delta_sum1       = xo.Float64[:]
    x_x_sum2         = xo.Float64[:]
    x_px_sum2        = xo.Float64[:]
    x_y_sum2         = xo.Float64[:]
    x_py_sum2        = xo.Float64[:]
    x_zeta_sum2      = xo.Float64[:]
    x_pzeta_sum2     = xo.Float64[:]
    x_delta_sum2     = xo.Float64[:]
    px_px_sum2       = xo.Float64[:]
    px_y_sum2        = xo.Float64[:]
    px_py_sum2       = xo.Float64[:]
    px_zeta_sum2     = xo.Float64[:]
    px_pzeta_sum2    = xo.Float64[:]
    px_delta_sum2    = xo.Float64[:]
    y_y_sum2         = xo.Float64[:]
    y_py_sum2        = xo.Float64[:]
    y_zeta_sum2      = xo.Float64[:]
    y_pzeta_sum2     = xo.Float64[:]
    y_delta_sum2     = xo.Float64[:]
    py_py_sum2       = xo.Float64[:]
    py_zeta_sum2     = xo.Float64[:]
    py_pzeta_sum2    = xo.Float64[:]
    py_delta_sum2    = xo.Float64[:]
    zeta_zeta_sum2   = xo.Float64[:]
    zeta_pzeta_sum2  = xo.Float64[:]
    zeta_delta_sum2  = xo.Float64[:]
    pzeta_pzeta_sum2 = xo.Float64[:]
    pzeta_delta_sum2 = xo.Float64[:]
    delta_delta_sum2 = xo.Float64[:]

_RECORD_FIELD_NAMES = {f.name for f in ParticleStatsMonitorRecord._fields}

class ParticleStatsMonitor(xt.BeamElement):
    _xofields = {
        'part_id_start':      xo.Int64,
        'part_id_end':        xo.Int64,
        'start_at_turn':      xo.Int64,
        'stop_at_turn':       xo.Int64,
        'frev':               xo.Float64,
        'sampling_frequency': xo.Float64,
        'cached':             xo.Int8[:],
        '_index':             xt.RecordIndex,
        'data':               ParticleStatsMonitorRecord,
        '_selector':          xo.Int16
    }

    _isthick = False
    needs_rng = False
    prototype = "monitor"
    allow_track = True
    iscollective = False
    has_backtrack = False
    behaves_like_drift = True
    allow_rot_and_shift = True
    allow_loss_refinement = False
    name_associated_aperture = None
    skip_in_loss_location_refinement = True

    _noexpr_fields = {'name', 'line'}
    _extra_c_sources = [
        '#include "xcoll/beam_elements/elements_src/monitor.h"'
    ]
    _store_in_to_dict = ["beta0", "gamma0", "mass0"]

    def __init__(self, **kwargs):
        """Monitor to save particle statistics (mean and variance of x,
        px, y, py, zeta, pzeta, delta) of (potentially a subset of)
        particles over a range of turns.

        The monitor allows for arbitrary sampling rate and can thus not
        only be used to monitor bunch positions, but also to record
        schottky spectra. Internally, the particle arrival time is used
        when determining the record index:

            i = sampling_frequency * (
                        ( at_turn - start_turn ) / f_rev
                        - zeta / beta0 / c0
                )

        where zeta=(s-beta0*c0*t) is the longitudinal coordinate of the
        particle, beta0 the relativistic beta factor of the particle, c0
        is the speed of light, at_turn is the current turn number, f_rev
        is the revolution frequency, and sampling_frequency is the
        sampling frequency.

        Note that the index is rounded, i.e. the result array represents
        data of particles equally distributed around the reference
        particle. For example, if the sampling_frequency is twice the
        revolution frequency, the first item contains data from
        particles in the range zeta/circumference = -0.25 .. 0.25, the
        second item in the range 0.25 .. 0.75 and so on.

        Args:
            num_particles (int, optional):
                Number of particles to monitor, starting from 0.
                Defaults to -1 which means ALL.
            particle_id_range (tuple, optional):
                Range of particle ids to monitor (start, stop). Stop is
                exclusive. Defaults to (particle_id_start,
                particle_id_start+num_particles).
            start_at_turn (int):
                First turn of reference particle (inclusive) at which to
                monitor. Defaults to 0.
            stop_at_turn (int):
                Last turn of reference particle (exclusive) at which to
                monitor. Defaults to start_at_turn + 1.
            frev (float):
                Revolution frequency in Hz of circulating beam (used to
                relate turn number to sample index). Defaults to 1.
            sampling_frequency (float):
                Sampling frequency in Hz. Defaults to 1.
            mean (bool):
                Whether or not to monitor means. Defaults to True.
            variance (bool):
                Whether or not to monitor variances. Defaults to True,
                unless mean is explicitly set to True.
            columns (list of str):
                List of variables to monitor. Can include any of 'x',
                'px', 'y', 'py', 'zeta', 'pzeta', 'delta'. Defaults to
                empty list, in which case all variables are monitored.
        """

        to_assign = {}
        if '_xobject' not in kwargs:
            # Get all flags related to what to monitor
            mon_mean = kwargs.pop('mean', None)
            mon_variance = kwargs.pop('variance', None)
            if mon_mean is None and mon_variance is None:
                mon_mean = mon_variance = True
            elif mon_mean is None:
                mon_mean = True
            elif mon_variance is None:
                mon_variance = not mon_mean
            columns = kwargs.pop('columns', [])
            mon_x = 'x' in columns
            mon_px = 'px' in columns
            mon_y = 'y' in columns
            mon_py = 'py' in columns
            mon_zeta = 'zeta' in columns
            mon_pzeta = 'pzeta' in columns
            mon_delta = 'delta' in columns
            if all(v is False for v in [mon_x, mon_px, mon_y, mon_py, mon_zeta,
                                        mon_pzeta, mon_delta]):
                mon_x = True
                mon_px = True
                mon_y = True
                mon_py = True
                mon_zeta = True
                mon_pzeta = True
                mon_delta = True
            kwargs['_selector']  =     int(mon_x)
            kwargs['_selector'] +=   2*int(mon_px)
            kwargs['_selector'] +=   4*int(mon_y)
            kwargs['_selector'] +=   8*int(mon_py)
            kwargs['_selector'] +=  16*int(mon_zeta)
            kwargs['_selector'] +=  32*int(mon_pzeta)
            kwargs['_selector'] +=  64*int(mon_delta)
            kwargs['_selector'] += 128*int(mon_mean)
            kwargs['_selector'] += 256*int(mon_variance)

            # Get ranges
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
                    kwargs['part_id_start'] = kwargs.pop(
                                                'particle_id_start', 0)
                    kwargs['part_id_end'] = kwargs['part_id_start'] \
                                            + num_particles
            start = int(kwargs.get('start_at_turn', 0))
            stop = int(kwargs.get('stop_at_turn', start + 1))
            frev = kwargs.get('frev', 1.)
            samp = kwargs.get('sampling_frequency', 1.)
            kwargs['start_at_turn'] = start
            kwargs['stop_at_turn'] = stop
            kwargs['frev'] = frev
            kwargs['sampling_frequency'] = samp

            # Prepare the arrays. Explicitely init with zeros (instead of size
            # only) to have consistent initial values
            if "data" not in kwargs:
                dd = {}
                size = int(round((stop-start) * samp / frev))
                for ff in ParticleStatsMonitorRecord._fields:
                    if ff.name == 'count':
                        dd[ff.name] = np.zeros(size, dtype=np.int64)
                    elif ff.name.endswith('_sum1'):
                        if not mon_mean:
                            dd[ff.name] = np.zeros(1, dtype=np.float64)
                        elif (ff.name.startswith('x_') and mon_x) or \
                             (ff.name.startswith('px_') and mon_px) or \
                             (ff.name.startswith('y_') and mon_y) or \
                             (ff.name.startswith('py_') and mon_py) or \
                             (ff.name.startswith('zeta_') and mon_zeta) or \
                             (ff.name.startswith('pzeta_') and mon_pzeta) or \
                             (ff.name.startswith('delta_') and mon_delta):
                            dd[ff.name] = np.zeros(size, dtype=np.float64)
                        else:
                            dd[ff.name] = np.zeros(1, dtype=np.float64)
                    elif ff.name.endswith('_sum2'):
                        if not mon_variance:
                            dd[ff.name] = np.zeros(1, dtype=np.float64)
                        else:
                            coords = ff.name[:-5].split('_')
                            if len(coords) != 2:
                                raise ValueError(f"Unexpected field {ff.name} "
                                              "in ParticleStatsMonitorRecord!")
                            if ('x' in coords and not mon_x) or \
                               ('px' in coords and not mon_px) or \
                               ('y' in coords and not mon_y) or \
                               ('py' in coords and not mon_py) or \
                               ('zeta' in coords and not mon_zeta) or \
                               ('pzeta' in coords and not mon_pzeta) or \
                               ('delta' in coords and not mon_delta):
                                dd[ff.name] = np.zeros(1, dtype=np.float64)
                            else:
                                dd[ff.name] = np.zeros(size, dtype=np.float64)
                    else:
                        raise ValueError(f"Unknown field {ff.name} in "
                                         "ParticleStatsMonitorRecord!")
                    kwargs['data'] = dd
            kwargs.setdefault('cached', np.zeros(len(kwargs['data']['count']),
                                                 dtype=np.int8))
            to_assign['_beta0'] = kwargs.pop('beta0', None)
            to_assign['_gamma0'] = kwargs.pop('gamma0', None)
            to_assign['_mass0'] = kwargs.pop('mass0', None)
        super().__init__(**kwargs)
        for kk, vv in to_assign.items():
            setattr(self, kk, vv)

    def to_json(self, file, indent=2, backend=['orjson', 'msgspec', 'json']):
        json.dump(self.to_dict(), file, indent=indent, backend=backend)

    @classmethod
    def from_json(cls, file):
        if not hasattr(file, '__iter__') or isinstance(file, (str, io.IOBase)):
            file = [file]
        dct = {}
        data = {}
        for f in file:
            this_dct = json.load(f)
            this_data = this_dct.pop('data')
            if dct == {}:
                dct = this_dct
                data = {kk: np.array(vv) for kk, vv in this_data.items()}
            else:
                if not deep_equal(this_dct, dct):
                    raise ValueError(f"JSON file {f} not compatible with "
                                      "previous ones!")
                for key, value in this_data.items():
                    if key not in data:
                        raise ValueError(f"JSON file {f} not compatible with "
                                          "previous ones!")
                    data[key] += np.array(value)
        self = cls.from_dict(dct | {'data': data})
        return self

    @classmethod
    def install(cls, line, name, *, at_s=None, at=None, s_tol=1.e-6, **kwargs):
        """Install a monitor in the line. The monitor will be configured
        with the line's reference particle parameters.
        """
        self = cls(**kwargs)
        if name in line.element_names:
            raise ValueError(f"Element {name} already exists in the line as "
                             f"{line[name].__class__.__name__}.")
        line.insert_element(element=self, name=name, at_s=at_s, at=at,
                            s_tol=s_tol)
        self.configure(line)
        return self

    def configure(self, line=None, *, beta0=None, gamma0=None, mass0=None):
        """Set optics parameters from line's reference particle. If
        `line` is not provided, `beta0`, `gamma0` and `mass0` must be
        provided explicitly.
        """
        if beta0 is not None and gamma0 is not None and mass0 is not None:
            self._beta0 = beta0
            self._gamma0 = gamma0
            self._mass0 = mass0
        elif line is not None and hasattr(line, 'particle_ref') \
        and line.particle_ref is not None:
            self._beta0 = line.particle_ref.beta0[0]
            self._gamma0 = line.particle_ref.gamma0[0]
            self._mass0 = line.particle_ref.mass0
        else:
            raise ValueError("Either a line with a particle_ref, or beta0, "
                             "gamma0 and mass0 must be provided!")

    def reset(self):
        """Reset the monitor data (to avoid unwanted accumulation)."""
        for field in [f.name for f in ParticleStatsMonitorRecord._fields]:
            ff = getattr(self.data, field)
            zeros = np.zeros(len(ff), dtype=ff._itemtype._dtype)
            setattr(self.data, field, zeros)
        for i in np.arange(len(self.count)):
            self.cached[i] = 0

    @property
    def beta0(self):
        return self._beta0

    @property
    def gamma0(self):
        return self._gamma0

    @property
    def mass0(self):
        return self._mass0

    @property
    def monitor_x(self):
        return bool(self._selector % 2)

    @property
    def monitor_px(self):
        return bool((self._selector >> 1) % 2)

    @property
    def monitor_y(self):
        return bool((self._selector >> 2) % 2)

    @property
    def monitor_py(self):
        return bool((self._selector >> 3) % 2)

    @property
    def monitor_zeta(self):
        return bool((self._selector >> 4) % 2)

    @property
    def monitor_pzeta(self):
        return bool((self._selector >> 5) % 2)

    @property
    def monitor_delta(self):
        return bool((self._selector >> 6) % 2)

    @property
    def monitor_mean(self):
        return bool((self._selector >> 7) % 2)

    @property
    def monitor_variance(self):
        return bool((self._selector >> 8) % 2)

    @property
    def turns(self):
        self._calculate()
        return self._turns

    @property
    def pc_mean(self):
        if self.monitor_delta and self.monitor_mean:
            self._calculate()
            one_plus_delta = 1 + self.delta_mean
            return one_plus_delta * self.beta0 * self.gamma0 * self.mass0
        else:
            raise ValueError("Momentum mean not available! Set "
                             "monitor_delta=True and monitor_mean=True.")

    @property
    def pc_var(self):
        if self.monitor_delta and self.monitor_variance:
            self._calculate()
            var_delta = self.delta_var
            return var_delta * self.beta0**2 * self.gamma0**2 * self.mass0**2
        else:
            raise ValueError("Momentum variance not available! Set "
                             "monitor_delta=True and monitor_variance=True.")

    @property
    def energy_mean(self):
        if self.monitor_pzeta and self.monitor_mean:
            self._calculate()
            one_plus_pzeta = 1 + self.beta0**2 * self.pzeta_mean
            return one_plus_pzeta * self.gamma0 * self.mass0
        else:
            raise ValueError("Energy mean not available! Set "
                             "monitor_pzeta=True and monitor_mean=True.")

    @property
    def energy_var(self):
        if self.monitor_pzeta and self.monitor_variance:
            self._calculate()
            var_pzeta = self.beta0**4 * self.pzeta_mean**2
            return var_pzeta * self.gamma0**2 * self.mass0**2
        else:
            raise ValueError("Energy variance not available! Set "
                             "monitor_pzeta=True and monitor_variance=True.")


    def _set_arr(self, name, mask, value):
        # Assign to attribute (create if does not exist yet):
        arr = getattr(self, name, None)
        if arr is None:
            arr = np.full(self.count.shape, np.nan, dtype=float)
            setattr(self, name, arr)
        arr[mask] = value

    def _calculate(self):
        N = self.count
        mask = (N > 0) & (self.cached == 0)
        if not np.any(mask):
            return

        # Calculate mean, variance, and std
        N = N[mask]
        self._turns = np.array(range(self.start_at_turn,
                                     self.stop_at_turn))[mask]

        # NaN for zero particles is expected behaviour
        with np.errstate(invalid='ignore'):
            for field in [f.name for f in ParticleStatsMonitorRecord._fields]:
                if field.endswith('_sum1'):
                    x = field[:-5]
                    ff = getattr(self, field)
                    ff = ff[mask] if len(ff) == len(mask) else ff
                    mean = ff / N
                    self._set_arr(f'_{x}_mean', mask, mean)
            for field in [f.name for f in ParticleStatsMonitorRecord._fields]:
                if field.endswith('_sum2'):
                    x1, x2 = field[:-5].split('_')
                    mean1 = getattr(self, f'_{x1}_mean')
                    mean2 = getattr(self, f'_{x2}_mean')
                    ff = getattr(self, field)
                    ff = ff[mask] if len(ff) == len(mask) else ff
                    variance = ff / (N - 1) - mean1 * mean2 * N / (N - 1)
                    self._set_arr(f'_{x1}_{x2}_var', mask, variance)
        for i in np.arange(len(self.count))[mask]:
            self.cached[i] = 1

    def __getattr__(self, attr):
        if attr in _RECORD_FIELD_NAMES:
            return getattr(self.data, attr).to_nparray()

        if attr.startswith('_'):
            raise AttributeError(f"Attribute {attr} not set!")

        if attr.endswith('_var') or attr.endswith('_mean'):
            self._calculate()
            return getattr(self, f'_{attr}')

        return super().__getattribute__(attr)


class EmittanceMonitor(ParticleStatsMonitor):
    _xofields = ParticleStatsMonitor._xofields
    _depends_on = [ParticleStatsMonitor]
    _extra_c_sources = [
        '#include "xcoll/beam_elements/elements_src/emittance_monitor.h"'
    ]

    def __init__(self, **kwargs):
        """
        Monitor to save the normalised beam emittance

        Similar to the BeamSizeMonitor and BeamPositionMonitor, it
        allows for arbitrary sampling rate and can thus not only be used
        to monitor bunch emittance, but also to record coasting beams.
        See their documentation for more information on how to use
        `frev` and `sampling_frequency`.

        Args:
            num_particles (int, optional): Number of particles to
                monitor, starting from 0. Defaults to -1 which means
                ALL.
            particle_id_range (tuple, optional): Range of particle ids
                to monitor (start, stop). Stop is exclusive. Defaults to
                (particle_id_start, particle_id_start+num_particles).
            start_at_turn (int): First turn of reference particle
                (inclusive) at which to monitor. Defaults to 0.
            stop_at_turn (int): Last turn of reference particle
                (exclusive) at which to monitor. Defaults to
                start_at_turn + 1.
            frev (float): Revolution frequency in Hz of circulating beam
                (used to relate turn number to sample index). Defaults
                to 1.
            sampling_frequency (float): Sampling frequency in Hz.
                Defaults to 1.
            horizontal (bool): Whether or not to monitor the horizontal
                plane. Defaults to True.
            vertical (bool): Whether or not to monitor the vertical
                plane. Defaults to True.
            longitudinal (bool): Whether or not to monitor the
                longitudinal plane. Defaults to True.
        """

        to_assign = {}
        if '_xobject' not in kwargs:
            kwargs['mean'] = True
            kwargs['variance'] = True
            horizontal = kwargs.pop('horizontal', None)
            vertical = kwargs.pop('vertical', None)
            longitudinal = kwargs.pop('longitudinal', None)
            # Small hack to allow e.g. horizontal=False to imply
            # vertical=True and longitudinal=True
            setval = np.unique([v for v in [horizontal, vertical, longitudinal]
                                if v is not None])
            if len(setval) == 0:
                setval = True
            elif len(setval) == 1:
                setval = not setval[0]
            else:
                setval = False
            if horizontal is None: horizontal = setval
            if vertical is None: vertical = setval
            if longitudinal is None: longitudinal = setval
            kwargs['columns'] = []
            if horizontal:
                kwargs['columns'] += ['x', 'px']
            if vertical:
                kwargs['columns'] += ['y', 'py']
            if longitudinal:
                kwargs['columns'] += ['zeta', 'pzeta']
            if kwargs.pop('monitor_delta', False):
                kwargs['columns'] += ['delta']
            to_assign['_suppress_warnings'] = kwargs.pop('suppress_warnings',
                                                         False)
        super().__init__(**kwargs)
        for kk, vv in to_assign.items():
            setattr(self, kk, vv)
        if not hasattr(self, 'cached_modes'):
            self.cached_modes = np.zeros(len(self.cached), dtype=np.int8)

    def reset(self):
        super().reset()
        self.cached_modes[:] = 0

    @property
    def suppress_warnings(self):
        return self._suppress_warnings

    @suppress_warnings.setter
    def suppress_warnings(self, val):
        self._suppress_warnings = bool(val)

    @property
    def gemitt_x(self):
        self._check_horizontal()
        self._calculate()
        return self._gemitt_x

    @property
    def gemitt_y(self):
        self._check_vertical()
        self._calculate()
        return self._gemitt_y

    @property
    def gemitt_zeta(self):
        self._check_longitudinal()
        self._calculate()
        return self._gemitt_zeta

    @property
    def nemitt_x(self):
        self._check_horizontal()
        self._calculate()
        return self._gemitt_x * self.beta0 * self.gamma0

    @property
    def nemitt_y(self):
        self._check_vertical()
        self._calculate()
        return self._gemitt_y * self.beta0 * self.gamma0

    @property
    def nemitt_zeta(self):
        self._check_longitudinal()
        self._calculate()
        return self._gemitt_zeta * self.beta0 * self.gamma0

    @property
    def gemitt_I(self):
        self._check_horizontal()
        self._calculate()
        self._calculate_modes()
        return self._gemitt_I

    @property
    def gemitt_II(self):
        self._check_vertical()
        self._calculate()
        self._calculate_modes()
        return self._gemitt_II

    @property
    def gemitt_III(self):
        self._check_longitudinal()
        self._calculate()
        self._calculate_modes()
        return self._gemitt_III

    @property
    def nemitt_I(self):
        self._check_horizontal()
        self._calculate()
        self._calculate_modes()
        return self._gemitt_I * self.beta0 * self.gamma0

    @property
    def nemitt_II(self):
        self._check_vertical()
        self._calculate()
        self._calculate_modes()
        return self._gemitt_II * self.beta0 * self.gamma0

    @property
    def nemitt_III(self):
        self._check_longitudinal()
        self._calculate()
        self._calculate_modes()
        return self._gemitt_III * self.beta0 * self.gamma0


    @property
    def horizontal(self):
        return self.monitor_x and self.monitor_px

    @property
    def vertical(self):
        return self.monitor_y and self.monitor_py

    @property
    def longitudinal(self):
        return self.monitor_zeta and self.monitor_pzeta


    def _check_horizontal(self):
        if not self.horizontal:
            raise ValueError("Horizontal plane not monitored!")

    def _check_vertical(self):
        if not self.vertical:
            raise ValueError("Vertical plane not monitored!")

    def _check_longitudinal(self):
        if not self.longitudinal:
            raise ValueError("Longitudinal plane not monitored!")


    def _calculate(self):
        N = self.count
        mask = (N > 0) & (self.cached == 0)
        if not np.any(mask):
            return

        super()._calculate()
        self.cached_modes[:] = 0

        # Calculate emittances
        gemitt_x = np.sqrt(self.x_x_var * self.px_px_var - self.x_px_var**2)
        gemitt_y = np.sqrt(self.y_y_var * self.py_py_var - self.y_py_var**2)
        gemitt_zeta = np.sqrt(self.zeta_zeta_var * self.pzeta_pzeta_var
                              - self.zeta_pzeta_var**2)
        setattr(self, '_gemitt_x', gemitt_x)
        setattr(self, '_gemitt_y', gemitt_y)
        setattr(self, '_gemitt_zeta', gemitt_zeta)


    def _calculate_modes(self):
        # Calculate emittance modes
        N = self.count
        mask = (N > 0) & (self.cached_modes == 0)
        if not np.any(mask):
            return

        S = np.array([[ 0., 1., 0., 0., 0., 0.],
                      [-1., 0., 0., 0., 0., 0.],
                      [ 0., 0., 0., 1., 0., 0.],
                      [ 0., 0.,-1., 0., 0., 0.],
                      [ 0., 0., 0., 0., 0., 1.],
                      [ 0., 0., 0., 0.,-1., 0.]])
        gemitt_I   = []
        gemitt_II  = []
        gemitt_III = []

        N = N[N > 0]
        for i in range(len(N)):
            if N[i] < 25:
                # Not enough statistics for a reliable calculation of the modes
                gemitt_I.append(0)
                gemitt_II.append(0)
                gemitt_III.append(0)
                continue

            if self.horizontal:
                block_x = np.array([[self.x_x_var[i],  self.x_px_var[i]],
                                    [self.x_px_var[i], self.px_px_var[i]]])
            else:
                block_x = np.zeros((2, 2))
            if self.vertical:
                block_y = np.array([[self.y_y_var[i],  self.y_py_var[i]],
                                    [self.y_py_var[i], self.py_py_var[i]]])
            else:
                block_y = np.zeros((2, 2))
            if self.longitudinal:
                block_z = np.array([[self.zeta_zeta_var[i],
                                     self.zeta_pzeta_var[i]],
                                    [self.zeta_pzeta_var[i],
                                     self.pzeta_pzeta_var[i]]])
            else:
                block_z = np.zeros((2, 2))
            if self.horizontal and self.vertical:
                block_xy = np.array([[self.x_y_var[i],  self.x_py_var[i]],
                                     [self.px_y_var[i], self.px_py_var[i]]])
            else:
                block_xy = np.zeros((2, 2))
            if self.horizontal and self.longitudinal:
                block_xz = np.array([[self.x_zeta_var[i],
                                      self.x_pzeta_var[i]],
                                     [self.px_zeta_var[i],
                                      self.px_pzeta_var[i]]])
            else:
                block_xz = np.zeros((2, 2))
            if self.vertical and self.longitudinal:
                block_yz = np.array([[self.y_zeta_var[i],
                                      self.y_pzeta_var[i]],
                                     [self.py_zeta_var[i],
                                      self.py_pzeta_var[i]]])
            else:
                block_yz = np.zeros((2, 2))

            covariance_S = np.dot(
                        np.block([[block_x,    block_xy,   block_xz],
                                  [block_xy.T, block_y,    block_yz],
                                  [block_xz.T, block_yz.T, block_z]]),
                        S)

            # Check for all zero matrix -> zero emittance
            if np.all(covariance_S < 1E-16):
                gemitt_I.append(0)
                gemitt_II.append(0)
                gemitt_III.append(0)
                continue

            cond_number = np.linalg.cond(covariance_S)
            if cond_number > 1e10 and not self.suppress_warnings:
                print(f"Warning: High condition number at time step {i}: "
                    + f"{cond_number}.\n{N[i]} particles logged.")

            rank = np.linalg.matrix_rank(covariance_S)
            expected_rank = int(self.horizontal) + int(self.vertical) + int(self.longitudinal)
            if rank < expected_rank and not self.suppress_warnings:
                print(f"Warning: Matrix is rank deficient at time step {i}: "
                    + f"rank {rank} instead of expected {len(covariance_S)}.\n"
                    + f"{N[i]} particles logged.")

            from xtrack.linear_normal_form import compute_linear_normal_form
            _, _, _, eigenvalues = compute_linear_normal_form(covariance_S)
            gemitt_I.append(eigenvalues[0].imag)
            gemitt_II.append(eigenvalues[1].imag)
            gemitt_III.append(eigenvalues[2].imag)

        setattr(self, '_gemitt_I',   np.array(gemitt_I))
        setattr(self, '_gemitt_II',  np.array(gemitt_II))
        setattr(self, '_gemitt_III', np.array(gemitt_III))
        self.cached_modes[:] = 1
