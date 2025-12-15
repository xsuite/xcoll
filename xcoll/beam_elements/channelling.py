# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xobjects as xo
import numpy as np
from ..general import _pkg_root
from .base import InvalidXcoll


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
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
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


class BentChannellingDev(InvalidXcoll):
    """
    Bent channelling element with ?harmonic-period step selection (needs to be done) and
    symplectic integrators.
    M2 = harmonic + nonlinear correction
    M3 = extended harmonic + bending correction
    M4 = Simplified Moliere + bending
    """

    _xofields = {

        # geometry / basic parameters
        'length' : xo.Float64,
        'U0'     : xo.Float64,     # potential depth (harmonic model)
        'Umax'   : xo.Float64,    # potential depth [eV]
        'R'      : xo.Float64,     # bending radius

        'dp'      : xo.Float64,    # interplanar distance [m]
        'aTF'     : xo.Float64,    # Thomas–Fermi screening length [m]
        'uT'      : xo.Float64,    # thermal vibration amplitude [m]

        'alpha_i' : xo.Float64,    # dimensionless
        'beta_i'  : xo.Float64,    # dimensionless

        # M2, M3, M4 are PHYSICAL MODELS
        'method'  : xo.Int8,       # 2, 3, 4

        # integrator configuration
        'order'   : xo.Int8,       # 2,4,6,8,10,12
        'variant' : xo.Int8,       # 1: Drift-Kick-Drift, 2: Kick-Drift-Kick

        # If <0 --> automatic step selection (Μ2, Μ3 --> 60/harmonic period and M4 --> 10 per harmonic period)
        'n_steps' : xo.Int64,
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements', 'elements_src', 'elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements', 'elements_src', 'bent_channelling_default_config.h'),
        _pkg_root.joinpath('beam_elements', 'elements_src', 'bent_channelling_kernels.h'),
        _pkg_root.joinpath('beam_elements', 'elements_src', 'bent_channelling_integrators.h'),
        _pkg_root.joinpath('beam_elements', 'elements_src', 'track_bent_channelling.h'),
        _pkg_root.joinpath('beam_elements', 'elements_src', 'bent_channelling.h'),


    ]

    # ------------------------------------------------------------------
    # Constructor with defaults (collimator-style)
    # ------------------------------------------------------------------
    

    def __init__(self, **kwargs):

        # =========================================================
        # Geometry
        # =========================================================
        kwargs.setdefault('length', 1e-4)
        kwargs.setdefault('R', 10.0)

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
        if 'n_steps' not in kwargs:
            method = kwargs['method']
            if method ==4:
                npp = 5
            else: 
                npp = 60 
            length = kwargs['length']
            Umax   = kwargs['Umax']
            alpha_i = kwargs['alpha_i']
            beta_i  = kwargs['beta_i']
            aTF     = kwargs['aTF']
            dp      = kwargs['dp']
            bpc     = kwargs.get('bpc', 150e9)  # safeguard

            Uxx0 = (
                2.0*Umax * alpha_i/beta_i
                * np.exp(-beta_i/aTF*0.5*dp)
                / bpc
                * (beta_i/aTF)**2
            )

            omega = np.sqrt(Uxx0)  # 1/m

            kwargs['n_steps'] = max(
                1,
                int(np.ceil(length * npp * omega / (2.0 * np.pi)))
            )

        # =========================================================
        # constructor LAST
        # =========================================================
        super().__init__(**kwargs)



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
        self.U0mol   = U0mol
        self.alpha_i = alpha_i
        self.beta_i  = beta_i
       
        
        
# Model 2
# --- Order 2 ---
class BentChannellingDevM2V1o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o02.h')
    ]
    
class BentChannellingDevM2V2o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o02.h')
    ] 
    
# --- Order 4 ---
class BentChannellingDevM2V1o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o04.h')
    ]
    
class BentChannellingDevM2V2o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o04.h')
    ] 
    
    
# --- Order 6 ---
class BentChannellingDevM2V1o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o06.h')
    ]
    
class BentChannellingDevM2V2o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o06.h')
    ] 
    
# --- Order 8 ---
class BentChannellingDevM2V1o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o08.h')
    ]
    
class BentChannellingDevM2V2o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o08.h')
    ]

# --- Order 10 ---
class BentChannellingDevM2V1o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o10.h')
    ]
    
class BentChannellingDevM2V2o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o10.h')
    ]

# --- Order 12 ---
class BentChannellingDevM2V1o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o12.h')
    ]
    
class BentChannellingDevM2V2o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o12.h')
    ]


# Model 3
# --- Order 2 ---
class BentChannellingDevM3V1o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o02.h')
    ]
    
class BentChannellingDevM3V2o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o02.h')
    ] 
    
# --- Order 4 ---
class BentChannellingDevM3V1o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o04.h')
    ]
    
class BentChannellingDevM3V2o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o04.h')
    ] 
    
    
# --- Order 6 ---
class BentChannellingDevM3V1o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o06.h')
    ]
    
class BentChannellingDevM3V2o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o06.h')
    ] 
    
# --- Order 8 ---
class BentChannellingDevM3V1o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o08.h')
    ]
    
class BentChannellingDevM3V2o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o08.h')
    ]

# --- Order 10 ---
class BentChannellingDevM3V1o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o10.h')
    ]
    
class BentChannellingDevM3V2o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o10.h')
    ]

# --- Order 12 ---
class BentChannellingDevM3V1o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o12.h')
    ]
    
class BentChannellingDevM3V2o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o12.h')
    ]

# Model 4
# --- Order 2 ---
class BentChannellingDevM4V1o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o02.h')
    ]
    
class BentChannellingDevM4V2o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o02.h')
    ] 
    
# --- Order 4 ---
class BentChannellingDevM4V1o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o04.h')
    ]
    
class BentChannellingDevM4V2o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o04.h')
    ] 
    
    
# --- Order 6 ---
class BentChannellingDevM4V1o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o06.h')
    ]
    
class BentChannellingDevM4V2o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o06.h')
    ] 
    
# --- Order 8 ---
class BentChannellingDevM4V1o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o08.h')
    ]
    
class BentChannellingDevM4V2o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o08.h')
    ]

# --- Order 10 ---
class BentChannellingDevM4V1o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o10.h')
    ]
    
class BentChannellingDevM4V2o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o10.h')
    ]

# --- Order 12 ---
class BentChannellingDevM4V1o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o12.h')
    ]
    
class BentChannellingDevM4V2o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
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
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o12.h')
    ]

                                       
