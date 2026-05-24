# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xobjects as xo
import xtrack as xt

from .base import BaseBlock, BaseCollimator, BaseCrystal, InvalidXcoll
from ..scattering_routines.everest import EverestEngine
from ..materials import Material, _DEFAULT_MATERIAL, _resolve_material
from ..general import _pkg_root


# TODO:
#      We want these elements to behave as if 'iscollective = True' when doing twiss etc (because they would ruin the CO),
#      but as if 'iscollective = False' for normal tracking as it is natively in C...
#      Currently this is achieved with the hack '_tracking' which defaults to False after installation in the line, and is
#      only activated around the track command. Furthermore, because of 'iscollective = False' we need to specify
#      get_backtrack_element. We want it nicer..


class EverestBlock(BaseBlock):
    _xofields = {**BaseBlock._xofields,
        '_material':        Material,
        'rutherford_rng':   xt.RandomRutherford,
        '_tracking':        xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    allow_double_sided = False
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = {*BaseBlock._noexpr_fields, 'material'}
    _skip_in_to_dict       = [*BaseBlock._skip_in_to_dict, '_material', '_tracking']
    _store_in_to_dict      = [*BaseBlock._store_in_to_dict, 'material']
    _internal_record_class = BaseBlock._internal_record_class

    _depends_on = [BaseBlock, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','everest_block.h')
    ]

    _kernels = {
        'EverestBlock_set_material': xo.Kernel(
                c_name='EverestBlock_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = _DEFAULT_MATERIAL
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


    @property
    def material(self):
        if self._material != _DEFAULT_MATERIAL:
            return self._material

    @material.setter
    def material(self, material):
        material = _resolve_material(material)
        if self.material != material:
            self._material = material
            self.EverestBlock_set_material(el=self)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                                 _buffer=_buffer, _offset=_offset)


class EverestCollimator(BaseCollimator):
    _xofields = {**BaseCollimator._xofields,
        '_material':        Material,
        'rutherford_rng':   xt.RandomRutherford,
        '_tracking':        xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    allow_double_sided = True
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = {*BaseCollimator._noexpr_fields, 'material'}
    _skip_in_to_dict       = [*BaseCollimator._skip_in_to_dict, '_material']
    _store_in_to_dict      = [*BaseCollimator._store_in_to_dict, 'material']
    _internal_record_class = BaseCollimator._internal_record_class

    _depends_on = [BaseCollimator, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','everest_collimator.h')
    ]

    _kernels = {
        'EverestCollimator_set_material': xo.Kernel(
                c_name='EverestCollimator_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = _DEFAULT_MATERIAL
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)

    @property
    def material(self):
        if self._material != _DEFAULT_MATERIAL:
            return self._material

    @material.setter
    def material(self, material):
        material = _resolve_material(material)
        if self.material != material:
            self._material = material
            self.EverestCollimator_set_material(el=self)

    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                            _buffer=_buffer, _offset=_offset)


class EverestCrystal(BaseCrystal):
    _xofields = {**BaseCrystal._xofields,
        'miscut':             xo.Float64,
        '_orient':            xo.Int8,
        '_critical_angle':    xo.Float64,
        '_critical_radius':   xo.Float64,
        '_material':          Material,
        'rutherford_rng':     xt.RandomRutherford,
        '_tracking':          xo.Int8
    }

    isthick = True
    needs_rng = True
    allow_track = True
    allow_double_sided = False
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = {*BaseCrystal._noexpr_fields, 'material', 'lattice'}
    _skip_in_to_dict       = [*BaseCrystal._skip_in_to_dict, '_orient', '_material']
    _store_in_to_dict      = [*BaseCrystal._store_in_to_dict, 'lattice', 'material']
    _internal_record_class = BaseCrystal._internal_record_class

    _depends_on = [BaseCrystal, EverestEngine]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','everest_crystal.h')
    ]

    _kernels = {
        'EverestCrystal_set_material': xo.Kernel(
                c_name='EverestCrystal_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            to_assign['material'] = kwargs.pop('material', None)
            kwargs['_material'] = _DEFAULT_MATERIAL
            to_assign['lattice'] = kwargs.pop('lattice', 'strip')
            kwargs.setdefault('miscut', 0)
            kwargs.setdefault('rutherford_rng', xt.RandomRutherford())
            kwargs.setdefault('_tracking', True)
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


    @property
    def material(self):
        if self._material != _DEFAULT_MATERIAL:
            return self._material

    @material.setter
    def material(self, material):
        material = _resolve_material(material, everest_crystal=True)
        if self.material != material:
            self._material = material
            self.EverestCrystal_set_material(el=self)

    @property
    def critical_angle(self):
        return self._critical_angle if abs(self._critical_angle) > 1.e-12 else None

    @property
    def critical_radius(self):
        return self._critical_radius if abs(self._critical_radius) > 1.e-10 else None

    @property
    def lattice(self):
        if self._orient == 1:
            return 'strip'
        elif self._orient == 2:
            return 'quasi-mosaic'
        else:
            raise ValueError(f"Illegal value {self._orient} for '_orient'!")

    @lattice.setter
    def lattice(self, lattice):
        if lattice == 'strip' or lattice == '110' or lattice == 110:
            self._orient = 1
        elif lattice == 'quasi-mosaic' or lattice == '111' or lattice == 111:
            self._orient = 2
        else:
            raise ValueError(f"Illegal value {lattice} for 'lattice'! "
                            + "Only use 'strip' (110) or 'quasi-mosaic' (111).")


    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidXcoll(length=-self.length, _context=_context,
                            _buffer=_buffer, _offset=_offset)


class EverestPreciseCrystal(EverestCrystal):
    """
    Bent channelling element with ?harmonic-period step selection (needs to be done) and
    symplectic integrators.
    M2 = harmonic + nonlinear correction
    M3 = extended harmonic + bending correction
    M4 = Simplified Moliere + bending
    """

    _xofields = {**EverestCrystal._xofields,
	    # material / potential parameters
        'aTF'     : xo.Float64,    # Thomas–Fermi screening length [m]
        'Umax'    : xo.Float64,     # potential depth [eV]
        'alpha_i' : xo.Float64,    # dimensionless
        'beta_i'  : xo.Float64,    # dimensionless

        # M2, M3, M4 are PHYSICAL MODELS
        'method'  : xo.Int8,       # 2, 3, 4

        # integrator configuration
        'order'   : xo.Int8,       # 2,4,6,8,10,12
        'variant' : xo.Int8,       # 1: Drift-Kick-Drift, 2: Kick-Drift-Kick

        # If <0 --> automatic step selection (Μ2, Μ3 --> 60/harmonic period and M4 --> 10 per harmonic period)
        'n_steps' : xo.Int64,
        '_n_steps_auto' : xo.Float64,
    }

    isthick = True
    needs_rng = True
    allow_track = True
    allow_double_sided = False
    behaves_like_drift = True
    allow_rot_and_shift = False
    allow_loss_refinement = True
    skip_in_loss_location_refinement = True

    _noexpr_fields         = EverestCrystal._noexpr_fields
    _skip_in_to_dict       = EverestCrystal._skip_in_to_dict
    _store_in_to_dict      = EverestCrystal._store_in_to_dict
    _internal_record_class = EverestCrystal._internal_record_class

    _depends_on = [EverestCrystal]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements', 'elements_src', 'everest_precise_crystal.h')
    ]

    _kernels = {
        'EverestPreciseCrystal_set_material': xo.Kernel(
                c_name='EverestPreciseCrystal_set_material',
                args=[xo.Arg(xo.ThisClass, name='el')]
            )
        }


    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            # =========================================================
            # Material defaults: Silicon (110)
            # =========================================================
            kwargs.setdefault('Umax',    23.9037)     # eV
            kwargs.setdefault('aTF',     1.94e-11)    # m
            kwargs.setdefault('alpha_i', 0.722452)
            kwargs.setdefault('beta_i',  0.573481)

            # =========================================================
            # Integrator defaults
            # =========================================================
            kwargs.setdefault('method',  4)   # M4
            kwargs.setdefault('order',   4)
            kwargs.setdefault('variant', 1)
            kwargs.setdefault('n_steps', -1)  # trigger for automatic step selection
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)
        self.update_auto_n_steps()


    @property
    def material(self):
        if self._material != _DEFAULT_MATERIAL:
            return self._material

    @material.setter
    def material(self, material):
        material = _resolve_material(material, everest_crystal=True)
        if self.material != material:
            self._material = material
            self.EverestPreciseCrystal_set_material(el=self)
            self.update_auto_n_steps()

    # user helpers
    @classmethod
    def get_available_methods(cls):
        return ['M2', 'M3', 'M4']

    def method_to_int(self, name):
        return {'M2': 2, 'M3': 3, 'M4': 4}[name]

    def set_method(self, method='M4', order=4, variant=1, n_steps=-1):
        self.method  = self.method_to_int(method)
        self.order   = order
        self.variant = variant
        self.n_steps = n_steps
        self.update_auto_n_steps()

    def set_simple_moliere_params(self, Umax, aTF, alpha_i, beta_i):
        """
        """
        self.Umax    = Umax
        self.aTF     = aTF
        self.alpha_i = alpha_i
        self.beta_i  = beta_i
        self.update_auto_n_steps()

    def update_auto_n_steps(self):
        if self.material:
            if self.method == 4:
                npp = 10
            else:
                npp = 60
            dp = self.material.crystal_plane_distance * 0.002 # [m]
            period = np.pi * dp * np.sqrt(0.5 / self.Umax)   # needs a factor sqrt(betapc)
            self._n_steps_auto = self.length / period * npp
