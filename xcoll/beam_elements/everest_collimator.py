import xobjects as xo
import xpart as xp

from .base_collimator import BaseCollimator, InvalidCollimator
from ..scattering_routines.everest import Material, CrystalMaterial
from ..general import _pkg_root



# TODO: 
#    - remove dx, dy, offset, tilt, as this should only be in colldb (and here only the jaw positions)
#    - We want these elements to behave as if 'iscollective = True' when doing twiss etc (because they would ruin the CO),
#      but as if 'iscollective = False' for normal tracking as it is natively in C...
#      Currently this is achieved with the hack '_tracking' which defaults to False after installation in the line, and is
#      only activated around the track command. Furthermore, because of 'iscollective = False' we need to specify get_backtrack_element
#      We want it nicer..

class EverestCollimator(BaseCollimator):
    _xofields = { **BaseCollimator._xofields,
        'offset':     xo.Float64,
        'onesided':   xo.Int8,
        'tilt':       xo.Float64[:],  # TODO: how to limit this to length 2
        'material':   Material,
        '_tracking':  xo.Int8
    }

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    _extra_c_sources = [
        xp.general._pkg_root.joinpath('random_number_generator/rng_src/base_rng.h'),
        xp.general._pkg_root.joinpath('random_number_generator/rng_src/local_particle_rng.h'),
        _pkg_root.joinpath('scattering_routines','everest','exponential_integral_Ei.h'),
        _pkg_root.joinpath('scattering_routines','everest','random.h'),
        _pkg_root.joinpath('scattering_routines','everest','scatter_init.h'),
        _pkg_root.joinpath('scattering_routines','everest','jaw.h'),
        _pkg_root.joinpath('scattering_routines','everest','scatter.h'),
        _pkg_root.joinpath('beam_elements','collimators_src','everest_collimator.h')
    ]

    _depends_on = [BaseCollimator, InvalidCollimator]

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
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

    # TODO: In principle we are not allowed to backtrack through a collimator
    #       However, the loss refinement will fail if this function is not provided
    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidCollimator(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)



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
        'material':    CrystalMaterial,
        '_tracking':   xo.Int8
    }

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    _extra_c_sources = [
        xp.general._pkg_root.joinpath('random_number_generator/rng_src/base_rng.h'),
        xp.general._pkg_root.joinpath('random_number_generator/rng_src/local_particle_rng.h'),
        _pkg_root.joinpath('scattering_routines','everest','exponential_integral_Ei.h'),
        _pkg_root.joinpath('scattering_routines','everest','random.h'),
        _pkg_root.joinpath('scattering_routines','everest','crystal.h'),
        _pkg_root.joinpath('scattering_routines','everest','scatter_crystal.h'),
        _pkg_root.joinpath('beam_elements','collimators_src','everest_crystal.h')
    ]

    _depends_on = [BaseCollimator, InvalidCollimator]

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            kwargs.setdefault('_tracking', True)
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

    # TODO: In principle we are not allowed to backtrack through a collimator
    #       However, the loss refinement will fail if this function is not provided
    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return InvalidCollimator(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)


