# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xobjects as xo
import numpy as np
from ..general import _pkg_root
from .base import InvalidXcoll, BaseCrystal


class ChannellingDev(InvalidXcoll):
    _xofields = {
        'length': xo.Float64,
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('scattering_routines','everest','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','channelling.h')
    ]

'''
class BentChannellingDev(InvalidXcoll):
    _xofields = {
        'length': xo.Float64,
        'method' : xo.Int64,  # 2, 3, 4
        'variant': xo.Int64,  # 1 or 2
        'order'  : xo.Int64,  # 2,4,6,8,10,12
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channelling.h')
    ]
               
'''


# what would be the proper class once the code is ready? xt.BeamElement with particular settings? BaseCrystal? something different?


class BentChannellingDev(BaseCrystal):
    """
    Bent channelling element with ?harmonic-period step selection (needs to be done) and
    symplectic integrators.
    M2 = harmonic + nonlinear correction
    M3 = extended harmonic + bending correction
    M4 = Simplified Moliere + bending
    """

    _xofields = {**BaseCrystal._xofields,

        

	 # material / potential parameters
        'dp'      : xo.Float64,    # interplanar distance [m]
        'aTF'     : xo.Float64,    # Thomas–Fermi screening length [m]
        'uT'      : xo.Float64,    # thermal vibration amplitude [m]


	'U0'     : xo.Float64,     # potential depth (harmonic model)
        'Umax'   : xo.Float64,    # potential depth [eV]
        
        'alpha_i' : xo.Float64,    # dimensionless
        'beta_i'  : xo.Float64,    # dimensionless

        # M2, M3, M4 are PHYSICAL MODELS
        'method'  : xo.Int8,       # 2, 3, 4

        # integrator configuration
        'order'   : xo.Int8,       # 2,4,6,8,10,12
        'variant' : xo.Int8,       # 1: Drift-Kick-Drift, 2: Kick-Drift-Kick

        # If <0 --> automatic step selection (Μ2, Μ3 --> 60/harmonic period and M4 --> 10 per harmonic period)
        'n_steps' : xo.Int64,
        '_n_steps_auto': xo.Float64,
        

    }

    isthick = True
    needs_rng = False #temporarilly
    allow_track = True
    allow_double_sided = False #no double sided
    behaves_like_drift = True
    allow_rot_and_shift = False #its own mechanism for these things
    allow_loss_refinement = True #for now only for compatibility reasons-->after integrating absorption will be useful
    skip_in_loss_location_refinement = True
    
    _noexpr_fields         = BaseCrystal._noexpr_fields
    _skip_in_to_dict       = BaseCrystal._skip_in_to_dict
    _store_in_to_dict      = BaseCrystal._store_in_to_dict
    _internal_record_class = BaseCrystal._internal_record_class

    _depends_on = [BaseCrystal]

    _extra_c_sources = [
        _pkg_root.joinpath('scattering_routines','everest', 'elliptic_functions.h'),
        _pkg_root.joinpath('scattering_routines','everest', 'precise_channelling_kernels.h'),
        _pkg_root.joinpath('scattering_routines','everest', 'precise_channelling_integrators.h'),
        _pkg_root.joinpath('beam_elements', 'elements_src', 'track_bent_channelling.h'),
        _pkg_root.joinpath('beam_elements', 'elements_src', 'bent_channelling.h'),


    ]

    # ------------------------------------------------------------------
    # Constructor with defaults (collimator-style)
    # ------------------------------------------------------------------
    

    def __init__(self, **kwargs):

        to_assign = {}

        if '_xobject' not in kwargs:
 
            # ---------------------------------------------------------
            # User-facing BaseCrystal geometry parameters
            # ---------------------------------------------------------
            required = ['length', 'width', 'height', 'bending_radius']

            missing = [key for key in required if key not in kwargs]
            if missing:
                raise ValueError(
                    "Missing required BentChannellingDev geometry parameters: "
                    + ", ".join(missing)
                    + ". Please provide length, width, height and bending_radius."
                )

            # Keep public parameters aside and assign them through properties
            # after BaseCrystal construction.
            to_assign['width'] = kwargs.pop('width')
            to_assign['height'] = kwargs.pop('height')
            to_assign['bending_radius'] = kwargs.pop('bending_radius')

            # Optional public parameters
            to_assign['angle'] = kwargs.pop('angle', 90.0)
            to_assign['tilt'] = kwargs.pop('tilt', 0.0)

            if 'gap' in kwargs:
                to_assign['gap'] = kwargs.pop('gap')

            if 'jaw_U' in kwargs:
                to_assign['jaw_U'] = kwargs.pop('jaw_U')

            # jaw_D is not stored directly in BaseCrystal, so ignore it for now.
            kwargs.pop('jaw_D', None)

            # BaseCrystal internal defaults only where necessary
            kwargs.setdefault('_side', 1)
            kwargs.setdefault('_align', 0)
            kwargs.setdefault('active', True)
            kwargs.setdefault('_record_interactions', False)
            kwargs.setdefault('_nemitt_x', 0.0)
            kwargs.setdefault('_nemitt_y', 0.0)
        
        # =========================================================
        # Material defaults: Silicon (110)
        # =========================================================
        kwargs.setdefault('U0',      21.7681)     # eV
        kwargs.setdefault('Umax',    23.9037)     # eV
        kwargs.setdefault('dp',      1.92e-10)    # m
        kwargs.setdefault('aTF',     1.94e-11)    # m
        kwargs.setdefault('uT',      7.5e-12)     # m
        kwargs.setdefault('alpha_i', 0.722452)
        kwargs.setdefault('beta_i',  0.573481)

        # =========================================================
        # Integrator defaults
        # =========================================================
        kwargs.setdefault('method',  4)   # M4
        kwargs.setdefault('order',   4)
        kwargs.setdefault('variant', 1)

        # how many steps per harmonic period
        #npp = kwargs.pop('n_steps_per_period', 60) I changed it to a method-sensitive one.

        # =========================================================
        # Automatic n_steps computation (ONLY if user did not set it)
        # =========================================================
        kwargs.setdefault('_n_steps_auto', 0.0)
        
        if 'n_steps' not in kwargs:
            kwargs['n_steps'] = -1  # trigger for automatic
            method = kwargs['method']
            if method ==4:
                npp = 10
            else: 
                npp = 60 
            length = kwargs['length']
            Umax   = kwargs['Umax']
            alpha_i = kwargs['alpha_i']
            beta_i  = kwargs['beta_i']
            aTF     = kwargs['aTF']
            dp      = kwargs['dp']
            # This value is only used to calculate the default number of steps.
            # The particle takes the real value from <track_bent_channelling.h>.
            #bpc     = kwargs.get('bpc', 150e9)  # safeguard

            Uxx0 = (
                2.0*Umax * alpha_i/beta_i
                * np.exp(-beta_i/aTF*0.5*dp)
            )

            omega = np.sqrt(Uxx0) * beta_i/aTF  # 1/m

            kwargs['_n_steps_auto'] = length * npp * omega / (2.0 * np.pi)
            
      
        # =========================================================
        # constructor LAST
        # =========================================================
        super().__init__(**kwargs)
        for key, val in to_assign.items():
            setattr(self, key, val)


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

    def set_simple_moliere_params(self,
                                  dp, aTF, uT,
                                  U0mol, alpha_i, beta_i):
        """
        Convenience helper to set the Simplified Molière potential parameters.
        Units:
          - dp, aTF, uT    in [m]
          - U0mol          in [eV]
          - alpha_i, beta_i dimensionless
        """
        self.dp      = dp
        self.aTF     = aTF
        self.uT      = uT
        self.U0mol   = U0mol #something's wrong here
        self.alpha_i = alpha_i
        self.beta_i  = beta_i
       
        
        

