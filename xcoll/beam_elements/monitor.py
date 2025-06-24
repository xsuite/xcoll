# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import io
import numpy as np

import xobjects as xo
import xtrack as xt

from ..general import _pkg_root


class ParticleStatsMonitorRecord(xo.Struct):
    count            = xo.Float64[:]
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
    px_px_sum2       = xo.Float64[:]
    px_y_sum2        = xo.Float64[:]
    px_py_sum2       = xo.Float64[:]
    px_zeta_sum2     = xo.Float64[:]
    px_pzeta_sum2    = xo.Float64[:]
    y_y_sum2         = xo.Float64[:]
    y_py_sum2        = xo.Float64[:]
    y_zeta_sum2      = xo.Float64[:]
    y_pzeta_sum2     = xo.Float64[:]
    py_py_sum2       = xo.Float64[:]
    py_zeta_sum2     = xo.Float64[:]
    py_pzeta_sum2    = xo.Float64[:]
    zeta_zeta_sum2   = xo.Float64[:]
    zeta_pzeta_sum2  = xo.Float64[:]
    pzeta_pzeta_sum2 = xo.Float64[:]
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
        '_plane_selector':    xo.Int8
    }

    behaves_like_drift = True
    allow_loss_refinement = True

    _noexpr_fields   = {'name', 'line'}
    _extra_c_sources = [
        xt._pkg_root.joinpath('headers/atomicadd.h'),
        _pkg_root.joinpath('beam_elements/elements_src/monitor.h')
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
            horizontal (bool): Whether or not to monitor the horizontal plane.
                Defaults to True.
            vertical (bool): Whether or not to monitor the vertical plane.
                Defaults to True.
            longitudinal (bool): Whether or not to monitor the longitudinal plane.
                Defaults to True.
            monitor_delta (bool): Whether or not to monitor delta and delta**2.
                Defaults to False.
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
            horizontal = kwargs.pop('horizontal', True)
            vertical = kwargs.pop('vertical', True)
            longitudinal = kwargs.pop('longitudinal', True)
            monitor_delta = kwargs.pop('monitor_delta', False)
            kwargs['_plane_selector']  =   int(horizontal)
            kwargs['_plane_selector'] += 2*int(vertical)
            kwargs['_plane_selector'] += 4*int(longitudinal)
            kwargs['_plane_selector'] += 8*int(monitor_delta)
            if "data" not in kwargs:
                # explicitely init with zeros (instead of size only) to have consistent initial values
                size = int(round((kwargs['stop_at_turn'] - kwargs['start_at_turn']) \
                                  * kwargs['sampling_frequency'] / kwargs['frev']))
                size_h = np.zeros(size) if horizontal else np.zeros(1)
                size_v = np.zeros(size) if vertical else np.zeros(1)
                size_l = np.zeros(size) if longitudinal else np.zeros(1)
                size_d = np.zeros(size) if monitor_delta else np.zeros(1)
                kwargs['data'] = {field.name: np.zeros(size) for field in ParticleStatsMonitorRecord._fields
                                  if 'x' not in field.name and 'y' not in field.name and 'zeta' not in field.name \
                                  and 'delta' not in field.name}
                kwargs['data'].update({field.name: size_h for field in ParticleStatsMonitorRecord._fields if 'x' in field.name})
                kwargs['data'].update({field.name: size_v for field in ParticleStatsMonitorRecord._fields if 'y' in field.name})
                kwargs['data'].update({field.name: size_l for field in ParticleStatsMonitorRecord._fields if 'zeta' in field.name})
                kwargs['data'].update({field.name: size_d for field in ParticleStatsMonitorRecord._fields if 'delta' in field.name})
        super().__init__(**kwargs)
        if not hasattr(self, '_cached'):
            self._cached = False
        if not hasattr(self, '_cached_modes'):
            self._cached_modes = False

    def to_json(self, file, indent=2):
        dct = self.to_dict()
        dct['data'] = self.data._to_json()
        dct['beta0'] = self.beta0
        dct['gamma0'] = self.gamma0
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
        self = cls.from_dict(dct | {'data': data})
        self._beta0 = beta0
        self._gamma0 = gamma0
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
            raise ValueError("Name not set! Install the monitor using xc.EmittanceMonitor.install() "
                             "or manually set the name after installation.")
        return self._name

    @name.setter
    def name(self, val):
        self._name = val

    @property
    def line(self):
        if not hasattr(self, '_line'):
            raise ValueError("Line not set! Install the monitor using xc.EmittanceMonitor.install() "
                             "or manually set the line after installation.")
        return self._line

    @line.setter
    def line(self, val):
        self._line = val


    @property
    def horizontal(self):
        return bool(self._plane_selector % 2)

    @property
    def vertical(self):
        return bool((self._plane_selector >> 1) % 2)

    @property
    def longitudinal(self):
        return bool((self._plane_selector >> 2) % 2)

    @property
    def monitor_delta(self):
        return bool((self._plane_selector >> 3) % 2)


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
    def turns(self):
        self._calculate()
        return self._turns


    def _calculate(self):
        if self._cached:
            return

        # Calculate mean and variance
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
        _pkg_root.joinpath('beam_elements/elements_src/emittance_monitor.h')
    ]

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
