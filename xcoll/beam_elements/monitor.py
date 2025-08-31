# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import io
import numpy as np

import xobjects as xo
import xtrack as xt


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

class ParticleStatsMonitor(xt.BeamElement):
    _xofields = {
        'part_id_start':      xo.Int64,
        'part_id_end':        xo.Int64,
        'start_at_turn':      xo.Int64,
        'stop_at_turn':       xo.Int64,
        'frev':               xo.Float64,
        'sampling_frequency': xo.Float64,
        '_index':             xt.RecordIndex,
        'data':               ParticleStatsMonitorRecord,
        '_cached':            xo.Int8,
        '_selector':          xo.Int16
    }

    behaves_like_drift = True
    allow_loss_refinement = True

    _noexpr_fields   = {'name', 'line'}
    _extra_c_sources = [
        '#include <xcoll/beam_elements/elements_src/monitor.h>'
    ]

    def __init__(self, **kwargs):
        """
        Monitor to save particle statistics (mean and variance of x, px, y, py,
        zeta, pzeta, delta) of (potentially a subset of) particles over a range
        of turns.

        The monitor allows for arbitrary sampling rate and can thus not only be
        used to monitor bunch positions, but also to record schottky spectra.
        Internally, the particle arrival time is used when determining the
        record index:

            i = sampling_frequency * ( ( at_turn - start_turn ) / f_rev - zeta / beta0 / c0 )

        where zeta=(s-beta0*c0*t) is the longitudinal coordinate of the
        particle, beta0 the relativistic beta factor of the particle, c0 is the
        speed of light, at_turn is the current turn number, f_rev is the
        revolution frequency, and sampling_frequency is the sampling frequency.

        Note that the index is rounded, i.e. the result array represents data
        of particles equally distributed around the reference particle. For
        example, if the sampling_frequency is twice the revolution frequency,
        the first item contains data from particles in the range
        zeta/circumference = -0.25 .. 0.25, the second item in the range
        0.25 .. 0.75 and so on.

        Args:
            num_particles (int, optional):
                Number of particles to monitor, starting from 0. Defaults to
                -1 which means ALL.
            particle_id_range (tuple, optional):
                Range of particle ids to monitor (start, stop). Stop is
                exclusive. Defaults to (particle_id_start, particle_id_start+
                num_particles).
            start_at_turn (int):
                First turn of reference particle (inclusive) at which to
                monitor. Defaults to 0.
            stop_at_turn (int):
                Last turn of reference particle (exclusive) at which to
                monitor. Defaults to start_at_turn + 1.
            frev (float):
                Revolution frequency in Hz of circulating beam (used to relate
                turn number to sample index). Defaults to 1.
            sampling_frequency (float):
                Sampling frequency in Hz. Defaults to 1.
            monitor_mean (bool):
                Whether or not to monitor means. Defaults to True.
            monitor_variance (bool):
                Whether or not to monitor variances. Defaults to True, unless
                monitor_mean is explicitly set to True.
            monitor_horizontal (bool):
                Whether or not to monitor x and px. Defaults to False, unless
                no coordinates are specified, in which case all are monitored.
            monitor_vertical (bool):
                Whether or not to monitor y and py. Defaults to False, unless
                no coordinates are specified, in which case all are monitored.
            monitor_longitudinal (bool):
                Whether or not to monitor zeta and pzeta. Defaults to False,
                unless no coordinates are specified, in which case all are
                monitored.
            monitor_x (bool):
                Whether or not to monitor x. Defaults to False, unless no
                coordinates are specified, in which case all are monitored.
            monitor_px (bool):
                Whether or not to monitor px. Defaults to False, unless no
                coordinates are specified, in which case all are monitored.
            monitor_y (bool):
                Whether or not to monitor y. Defaults to False, unless no
                coordinates are specified, in which case all are monitored.
            monitor_py (bool):
                Whether or not to monitor py. Defaults to False, unless no
                coordinates are specified, in which case all are monitored.
            monitor_zeta (bool):
                Whether or not to monitor zeta. Defaults to False, unless no
                coordinates are specified, in which case all are monitored.
            monitor_pzeta (bool):
                Whether or not to monitor pzeta. Defaults to False, unless no
                coordinates are specified, in which case all are monitored.
            monitor_delta (bool):
                Whether or not to monitor delta. Defaults to False, unless no
                coordinates are specified, in which case all are monitored.
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

            # Get all flags related to what to monitor
            monitor_x = kwargs.pop('monitor_x', None)
            monitor_px = kwargs.pop('monitor_px', None)
            monitor_y = kwargs.pop('monitor_y', None)
            monitor_py = kwargs.pop('monitor_py', None)
            monitor_zeta = kwargs.pop('monitor_zeta', None)
            monitor_pzeta = kwargs.pop('monitor_pzeta', None)
            monitor_delta = kwargs.pop('monitor_delta', None)
            monitor_horizontal = kwargs.pop('monitor_horizontal', None)
            monitor_vertical = kwargs.pop('monitor_vertical', None)
            monitor_longitudinal = kwargs.pop('monitor_longitudinal', None)
            if monitor_horizontal is not None:
                if monitor_x is not None:
                    raise ValueError("Cannot specify both monitor_horizontal and monitor_x!")
                if monitor_px is not None:
                    raise ValueError("Cannot specify both monitor_horizontal and monitor_px!")
                monitor_x = monitor_horizontal
                monitor_px = monitor_horizontal
            if monitor_vertical is not None:
                if monitor_y is not None:
                    raise ValueError("Cannot specify both monitor_vertical and monitor_y!")
                if monitor_py is not None:
                    raise ValueError("Cannot specify both monitor_vertical and monitor_py!")
                monitor_y = monitor_vertical
                monitor_py = monitor_vertical
            if monitor_longitudinal is not None:
                if monitor_zeta is not None:
                    raise ValueError("Cannot specify both monitor_longitudinal and monitor_zeta!")
                if monitor_pzeta is not None:
                    raise ValueError("Cannot specify both monitor_longitudinal and monitor_pzeta!")
                monitor_zeta = monitor_longitudinal
                monitor_pzeta = monitor_longitudinal
                monitor_delta = monitor_longitudinal
            if all(v is None for v in [monitor_x, monitor_px, monitor_y, monitor_py,
                                       monitor_zeta, monitor_pzeta, monitor_delta]):
                monitor_x = True
                monitor_px = True
                monitor_y = True
                monitor_py = True
                monitor_zeta = True
                monitor_pzeta = True
                monitor_delta = True
            else:
                # This little hack allows to do monitor_delta = False implying all else should be True
                setval = np.unique([v for v in [monitor_x, monitor_px, monitor_y, monitor_py,
                                               monitor_zeta, monitor_pzeta, monitor_delta]
                                   if v is not None])
                if len(setval) > 1:
                    setval = False
                else:
                    setval = not setval[0]
                if monitor_x is None: monitor_x = setval
                if monitor_px is None: monitor_px = setval
                if monitor_y is None: monitor_y = setval
                if monitor_py is None: monitor_py = setval
                if monitor_zeta is None: monitor_zeta = setval
                if monitor_pzeta is None: monitor_pzeta = setval
                if monitor_delta is None: monitor_delta = setval
            if all(v is False for v in [monitor_x, monitor_px, monitor_y, monitor_py,
                                        monitor_zeta, monitor_pzeta, monitor_delta]):
                raise ValueError("At least one of the coordinates must be monitored!")
            monitor_mean = kwargs.pop('monitor_mean', None)
            monitor_variance = kwargs.pop('monitor_variance', None)
            if monitor_mean is None and monitor_variance is None:
                monitor_mean = True
                monitor_variance = True
            elif monitor_mean is None:
                monitor_mean = True
            elif monitor_variance is None:
                monitor_variance = not monitor_mean
            kwargs['_selector']  =     int(monitor_x)
            kwargs['_selector'] +=   2*int(monitor_px)
            kwargs['_selector'] +=   4*int(monitor_y)
            kwargs['_selector'] +=   8*int(monitor_py)
            kwargs['_selector'] +=  16*int(monitor_zeta)
            kwargs['_selector'] +=  32*int(monitor_pzeta)
            kwargs['_selector'] +=  64*int(monitor_delta)
            kwargs['_selector'] += 128*int(monitor_mean)
            kwargs['_selector'] += 256*int(monitor_variance)

            # Prepare the arrays. Explicitely init with zeros (instead of size only) to have consistent initial values
            if "data" not in kwargs:
                size = int(round((kwargs['stop_at_turn'] - kwargs['start_at_turn']) \
                                  * kwargs['sampling_frequency'] / kwargs['frev']))
                kwargs['data'] = {}
                for field in ParticleStatsMonitorRecord._fields:
                    if field.name == 'count':
                        kwargs['data'].update({field.name: np.zeros(size, dtype=np.int64)})
                    elif field.name.endswith('_sum1'):
                        if not monitor_mean:
                            kwargs['data'].update({field.name: np.zeros(1, dtype=np.float64)})
                        elif (field.name.startswith('x_') and monitor_x) or \
                             (field.name.startswith('px_') and monitor_px) or \
                             (field.name.startswith('y_') and monitor_y) or \
                             (field.name.startswith('py_') and monitor_py) or \
                             (field.name.startswith('zeta_') and monitor_zeta) or \
                             (field.name.startswith('pzeta_') and monitor_pzeta) or \
                             (field.name.startswith('delta_') and monitor_delta):
                            kwargs['data'].update({field.name: np.zeros(size, dtype=np.float64)})
                        else:
                            kwargs['data'].update({field.name: np.zeros(1, dtype=np.float64)})
                    elif field.name.endswith('_sum2'):
                        if not monitor_variance:
                            kwargs['data'].update({field.name: np.zeros(1, dtype=np.float64)})
                        else:
                            coords = field.name[:-5].split('_')
                            if len(coords) != 2:
                                raise ValueError(f"Unexpected field {field.name} in ParticleStatsMonitorRecord!")
                            if ('x' in coords and not monitor_x) or \
                               ('px' in coords and not monitor_px) or \
                               ('y' in coords and not monitor_y) or \
                               ('py' in coords and not monitor_py) or \
                               ('zeta' in coords and not monitor_zeta) or \
                               ('pzeta' in coords and not monitor_pzeta) or \
                               ('delta' in coords and not monitor_delta):
                                kwargs['data'].update({field.name: np.zeros(1, dtype=np.float64)})
                            else:
                                kwargs['data'].update({field.name: np.zeros(size, dtype=np.float64)})
                    else:
                        raise ValueError(f"Unknown field {field.name} in ParticleStatsMonitorRecord!")
        super().__init__(**kwargs)
        if not hasattr(self, '_cached'):
            self._cached = False

    def to_json(self, file, indent=2):
        dct = self.to_dict()
        dct['data'] = self.data._to_json()
        dct['beta0'] = self.beta0
        dct['gamma0'] = self.gamma0
        dct['mass0'] = self.mass0
        xt.json.dump(dct, file, indent=indent)

    @classmethod
    def from_json(cls, file):
        if not hasattr(file, '__iter__') or isinstance(file, (str, io.IOBase)):
            file = [file]
        dct = {}
        data = {}
        for f in file:
            this_dct = xt.json.load(f)
            this_data = this_dct.pop('data')
            if dct == {}:
                dct = this_dct
                data = {kk: np.array(vv) for kk, vv in this_data.items()}
            else:
                if not xt.line._dicts_equal(this_dct, dct):
                    raise ValueError(f"Json file {f} not compatible with previous ones!")
                for key, value in this_data.items():
                    if key not in data:
                        raise ValueError(f"Json file {f} not compatible with previous ones!")
                    data[key] += np.array(value)
        beta0 = dct.pop('beta0')
        gamma0 = dct.pop('gamma0')
        mass0 = dct.pop('mass0')
        self = cls.from_dict(dct | {'data': data})
        self._beta0 = beta0
        self._gamma0 = gamma0
        self._mass0 = mass0
        return self

    @classmethod
    def install(cls, line, name, *, at_s=None, at=None, s_tol=1.e-6, **kwargs):
        self = cls(**kwargs)
        if name in line.element_names:
            raise ValueError(f"Element {name} already exists in the line as {line[name].__class__.__name__}.")
        line.insert_element(element=self, name=name, at_s=at_s, at=at, s_tol=s_tol)
        self._name = name
        self._line = line
        return self

    @property
    def name(self):
        if not hasattr(self, '_name'):
            raise ValueError(f"Name not set! Install the monitor using {self.__class__.__name__}.install() "
                              "or manually set the name after installation.")
        return self._name

    @name.setter
    def name(self, val):
        self._name = val

    @property
    def line(self):
        if not hasattr(self, '_line'):
            raise ValueError(f"Line not set! Install the monitor using {self.__class__.__name__}.install() "
                              "or manually set the line after installation.")
        return self._line

    @line.setter
    def line(self, val):
        self._line = val

    @property
    def beta0(self):
        if not hasattr(self, '_beta0'):
            self._beta0 = self.line.particle_ref.beta0[0]
        return self._beta0

    @property
    def gamma0(self):
        if not hasattr(self, '_gamma0'):
            self._gamma0 = self.line.particle_ref.gamma0[0]
        return self._gamma0

    @property
    def mass0(self):
        if not hasattr(self, '_mass0'):
            self._mass0 = self.line.particle_ref.mass0
        return self._mass0

    # TODO: need to store mass_ratio!

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
    def horizontal(self):
        return self.monitor_x and self.monitor_px

    @property
    def vertical(self):
        return self.monitor_y and self.monitor_py

    @property
    def longitudinal(self):
        return self.monitor_zeta and self.monitor_pzeta

    @property
    def turns(self):
        self._calculate()
        return self._turns

    @property
    def pc_mean(self):
        if self.monitor_delta and self.monitor_mean:
            self._calculate()
            return (1 + self.delta_mean) * self.beta0 * self.gamma0 * self.mass0
        else:
            raise ValueError("Momentum mean not available! Set monitor_delta=True and monitor_mean=True.")

    @property
    def pc_var(self):
        if self.monitor_delta and self.monitor_mean:
            self._calculate()
            return self.delta_var * self.beta0**2 * self.gamma0**2 * self.mass0**2
        else:
            raise ValueError("Momentum variance not available! Set monitor_delta=True and monitor_variance=True.")

    @property
    def energy_mean(self):
        if self.monitor_pzeta and self.monitor_mean:
            self._calculate()
            return (1 + self.beta0**2 * self.pzeta_mean) * self.gamma0 * self.mass0
        else:
            raise ValueError("Energy mean not available! Set monitor_pzeta=True and monitor_mean=True.")

    @property
    def energy_var(self):
        if self.monitor_pzeta and self.monitor_variance:
            self._calculate()
            return self.beta0**4 * self.pzeta_mean**2 * self.gamma0**2 * self.mass0**2
        else:
            raise ValueError("Energy variance not available! Set monitor_pzeta=True and monitor_variance=True.")


    def _calculate(self):
        if self._cached:
            return

        N = self.count
        mask = N > 0
        N = N[mask]
        self._turns = np.array(range(self.start_at_turn, self.stop_at_turn))[mask]
        with np.errstate(invalid='ignore'):  # NaN for zero particles is expected behaviour
            for field in [f.name for f in ParticleStatsMonitorRecord._fields]:
                if field.endswith('_sum1'):
                    x = field[:-5]
                    ff = getattr(self, f'{x}_sum1')
                    ff = ff[mask] if len(ff) == len(mask) else ff
                    mean = ff / N
                    setattr(self, f'_{x}_mean', mean)
            for field in [f.name for f in ParticleStatsMonitorRecord._fields]:
                if field.endswith('_sum2'):
                    x1, x2 = field[:-5].split('_')
                    mean1 = getattr(self, f'_{x1}_mean')
                    mean2 = getattr(self, f'_{x2}_mean')
                    ff = getattr(self, field)
                    ff = ff[mask] if len(ff) == len(mask) else ff
                    variance = ff / (N - 1) - mean1 * mean2 * N / (N - 1)
                    setattr(self, f'_{x1}_{x2}_var', variance)
        self._cached = True


    def __getattr__(self, attr):
        if attr in [f.name for f in ParticleStatsMonitorRecord._fields]:
            return getattr(self.data, attr).to_nparray()

        elif attr in self.__class__._xofields:
            return super().__getattr(attr)

        else:
            if attr.startswith('_'):
                raise AttributeError(f"Attribute {attr} not set!")

            if attr.endswith('_var'):
                self._calculate()
                return getattr(self, f'_{attr}')

            elif attr.endswith('_mean'):
                self._calculate()
                return getattr(self, f'_{attr}')

            else:
                raise AttributeError(f"{self.__class__.__name__} has no attribute '{attr}'")


class EmittanceMonitor(ParticleStatsMonitor):
    _xofields = ParticleStatsMonitor._xofields

    behaves_like_drift = True
    allow_loss_refinement = True

    _depends_on      = [ParticleStatsMonitor]
    _noexpr_fields   = ParticleStatsMonitor._noexpr_fields
    _extra_c_sources = [
        '#include <xcoll/beam_elements/elements_src/emittance_monitor.h>'
    ]

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs['monitor_mean'] = True
            kwargs['monitor_variance'] = True
            kwargs['monitor_horizontal'] = kwargs.pop('horizontal', None)
            kwargs['monitor_vertical'] = kwargs.pop('vertical', None)
            kwargs['monitor_longitudinal'] = kwargs.pop('longitudinal', None)
        super().__init__(**kwargs)
        if not hasattr(self, '_cached_modes'):
            self._cached_modes = False

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
        if self._cached:
            return

        super()._calculate()
        self._cached_modes = False

        # Calculate emittances
        gemitt_x = np.sqrt(self.x_x_var * self.px_px_var - self.x_px_var**2)
        gemitt_y = np.sqrt(self.y_y_var * self.py_py_var - self.y_py_var**2)
        gemitt_zeta = np.sqrt(self.zeta_zeta_var * self.pzeta_pzeta_var - self.zeta_pzeta_var**2)
        setattr(self, '_gemitt_x', gemitt_x)
        setattr(self, '_gemitt_y', gemitt_y)
        setattr(self, '_gemitt_zeta', gemitt_zeta)


    def _calculate_modes(self):
        if self._cached_modes:
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
        N = self.count
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
                block_z = np.array([[self.zeta_zeta_var[i],  self.zeta_pzeta_var[i]],
                                    [self.zeta_pzeta_var[i], self.pzeta_pzeta_var[i]]])
            else:
                block_z = np.zeros((2, 2))
            if self.horizontal and self.vertical:
                block_xy = np.array([[self.x_y_var[i],  self.x_py_var[i]],
                                     [self.px_y_var[i], self.px_py_var[i]]])
            else:
                block_xy = np.zeros((2, 2))
            if self.horizontal and self.longitudinal:
                block_xz = np.array([[self.x_zeta_var[i],  self.x_pzeta_var[i]],
                                     [self.px_zeta_var[i], self.px_pzeta_var[i]]])
            else:
                block_xz = np.zeros((2, 2))
            if self.vertical and self.longitudinal:
                block_yz = np.array([[self.y_zeta_var[i],  self.y_pzeta_var[i]],
                                     [self.py_zeta_var[i], self.py_pzeta_var[i]]])
            else:
                block_yz = np.zeros((2, 2))

            covariance_S = np.dot(np.block([[block_x,    block_xy,   block_xz],
                                            [block_xy.T, block_y,    block_yz],
                                            [block_xz.T, block_yz.T, block_z]]),
                                  S)
            cond_number = np.linalg.cond(covariance_S)
            if cond_number > 1e10:
                print(f"Warning: High condition number when calculating "
                    + f"the emittances modes at time step {i}: {cond_number}.\n"
                    + f"One of the coordinates might be close to zero or not "
                    + f"varying enough among the different particles. Only "
                    + f"{N[i]} particles were logged at this step.")
            rank = np.linalg.matrix_rank(covariance_S)
            expected_rank = int(self.horizontal) + int(self.vertical) + int(self.longitudinal)
            if rank < expected_rank:
                print(f"Warning: Matrix is rank deficient when calculating "
                    + f"the emittances modes at time step {i}: rank {rank} "
                    + f"instead of expected {len(covariance_S)}.\n"
                    + f"One of the coordinates might be close to zero or not "
                    + f"varying enough among the different particles. Only "
                    + f"{N[i]} particles were logged at this step.")

            from xtrack.linear_normal_form import compute_linear_normal_form
            _, _, _, eigenvalues = compute_linear_normal_form(covariance_S)
            gemitt_I.append(eigenvalues[0].imag)
            gemitt_II.append(eigenvalues[1].imag)
            gemitt_III.append(eigenvalues[2].imag)

        setattr(self, '_gemitt_I',   np.array(gemitt_I))
        setattr(self, '_gemitt_II',  np.array(gemitt_II))
        setattr(self, '_gemitt_III', np.array(gemitt_III))
        self._cached_modes = True
