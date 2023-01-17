import xobjects as xo

from .base_collimator import BaseCollimator
from ..scattering_routines.everest import Material, CrystalMaterial
from ..general import _pkg_root



# TODO: remove dx, dy, offset, tilt, as this should only be in colldb (and here only the jaw positions)
class EverestCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'offset':     xo.Float64,
        'onesided':   xo.Int8,
        'tilt':       xo.Float64[:],  # TODO: how to limit this to length 2
        'material':   Material
    }

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    iscollective = False

    _extra_c_sources = [
        _pkg_root.joinpath('scattering_routines','everest','exponential_integral_Ei.h'),
        _pkg_root.joinpath('scattering_routines','everest','random.h'),
        _pkg_root.joinpath('scattering_routines','everest','scatter_init.h'),
        _pkg_root.joinpath('scattering_routines','everest','jaw.h'),
        _pkg_root.joinpath('scattering_routines','everest','scatter.h'),
        _pkg_root.joinpath('beam_elements','collimators_src','everest_collimator.h')
    ]

    _depends_on = [BaseCollimator]

    def __init__(self, **kwargs):
        kwargs.setdefault('offset', 0)
        kwargs.setdefault('onesided', False)
        kwargs.setdefault('tilt', [0,0])
        tilt = kwargs['tilt']
        if hasattr(tilt, '__iter__'):
            if isinstance(tilt, str):
                raise ValueError("Variable tilt has to be a number or array of numbers!")
            elif len(tilt) == 1:
                tilt = [tilt[0], tilt[0]]
            elif len(tilt) > 2:
                raise ValueError("Variable tilt cannot have more than two elements (tilt_L and tilt_R)!")
        else:
            tilt = [tilt, tilt]
        kwargs['tilt'] = tilt
        super().__init__(**kwargs)



class EverestCrystal(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'align_angle': xo.Float64,  #  = - sqrt(eps/beta)*alpha*nsigma
        'bend':        xo.Float64,
        'xdim':        xo.Float64,
        'ydim':        xo.Float64,
        'thick':       xo.Float64,
        'crytilt':     xo.Float64,
        'miscut':      xo.Float64,
        'orient':      xo.Float64,
        'offset':      xo.Float64,
        'onesided':    xo.Int8,
        'tilt':        xo.Float64[:],  # TODO: how to limit this to length 2
        'material':    CrystalMaterial
    }

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    iscollective = False

    _extra_c_sources = [
        _pkg_root.joinpath('scattering_routines','everest','exponential_integral_Ei.h'),
        _pkg_root.joinpath('scattering_routines','everest','random.h'),
        _pkg_root.joinpath('scattering_routines','everest','crystal.h'),
        _pkg_root.joinpath('scattering_routines','everest','scatter_crystal.h'),
        _pkg_root.joinpath('beam_elements','collimators_src','everest_crystal.h')
    ]

    _depends_on = [BaseCollimator]

    def __init__(self, **kwargs):
        kwargs.setdefault('offset', 0)
        kwargs.setdefault('onesided', False)
        kwargs.setdefault('tilt', [0,0])
        tilt = kwargs['tilt']
        if hasattr(tilt, '__iter__'):
            if isinstance(tilt, str):
                raise ValueError("Variable tilt has to be a number or array of numbers!")
            elif len(tilt) == 1:
                tilt = [tilt[0], tilt[0]]
            elif len(tilt) > 2:
                raise ValueError("Variable tilt cannot have more than two elements (tilt_L and tilt_R)!")
        else:
            tilt = [tilt, tilt]
        kwargs['tilt'] = tilt
        kwargs.setdefault('bend', 0)
        kwargs.setdefault('xdim', 0)
        kwargs.setdefault('ydim', 0)
        kwargs.setdefault('thick', 0)
        kwargs.setdefault('crytilt', 0)
        kwargs.setdefault('miscut', 0)
        kwargs.setdefault('orient', 0)
        super().__init__(**kwargs)


