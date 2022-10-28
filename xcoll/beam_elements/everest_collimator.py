import numpy as np

import xobjects as xo

from .collimators import BaseCollimator
from ..scattering_routines.everest import track, Material, CrystalMaterial


# TODO: remove dx, dy, offset, tilt, as this should only be in colldb (and here only the jaw positions)
class Collimator(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        'dpx':        xo.Float64,
        'dpy':        xo.Float64,
        'offset':     xo.Float64,
        'onesided':   xo.Int8,
        'tilt':       xo.Float64[:],  # TODO: how to limit this to length 2
        'material':   Material
    }

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    iscollective = True # TODO: will be set to False when fully in C

    def __init__(self, **kwargs):
        kwargs.setdefault('dpx', 0)
        kwargs.setdefault('dpx', 0)
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


    def track(self, particles):  # TODO: write impacts
        npart = particles._num_active_particles
        if npart == 0:
            return
        if self.material is None:
            raise ValueError("Cannot track if material is not set!")
        
        # TODO: when in C, drifting should call routine from drift element
        #       such that can be exact etc
        if not self.is_active:
            # Drift full length
            L = self.length
            if L > 0:
                rpp = particles.rpp[:npart]
                xp = particles.px[:npart] * rpp
                yp = particles.py[:npart] * rpp
                dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
                particles.x[:npart] += xp * L
                particles.y[:npart] += yp * L
                particles.s[:npart] += L
                particles.zeta[:npart] += dzeta*L
        else:
            # Drift inactive front
            L = self.inactive_front
            if L > 0:
                rpp = particles.rpp[:npart]
                xp = particles.px[:npart] * rpp
                yp = particles.py[:npart] * rpp
                dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
                particles.x[:npart] += xp * L
                particles.y[:npart] += yp * L
                particles.s[:npart] += L
                particles.zeta[:npart] += dzeta*L
  
            track(self, particles, npart)

            # Drift inactive back
            L = self.inactive_back
            if L > 0:
                rpp = particles.rpp[:npart]
                xp = particles.px[:npart] * rpp
                yp = particles.py[:npart] * rpp
                dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
                particles.x[:npart] += xp * L
                particles.y[:npart] += yp * L
                particles.s[:npart] += L
                particles.zeta[:npart] += dzeta*L

            particles.reorganize()



class Crystal(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        'dpx':         xo.Float64,
        'dpy':         xo.Float64,
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

    iscollective = True # TODO: will be set to False when fully in C

    def __init__(self, **kwargs):
        kwargs.setdefault('dpx', 0)
        kwargs.setdefault('dpx', 0)
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


    def track(self, particles):  # TODO: write impacts
        npart = particles._num_active_particles
        if npart == 0:
            return
        if self.material is None:
            raise ValueError("Cannot track if material is not set!")
        
        # TODO: when in C, drifting should call routine from drift element
        #       such that can be exact etc
        if not self.is_active:
            # Drift full length
            L = self.length
            if L > 0:
                rpp = particles.rpp[:npart]
                xp = particles.px[:npart] * rpp
                yp = particles.py[:npart] * rpp
                dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
                particles.x[:npart] += xp * L
                particles.y[:npart] += yp * L
                particles.s[:npart] += L
                particles.zeta[:npart] += dzeta*L
        else:
            # Drift inactive front
            L = self.inactive_front
            if L > 0:
                rpp = particles.rpp[:npart]
                xp = particles.px[:npart] * rpp
                yp = particles.py[:npart] * rpp
                dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
                particles.x[:npart] += xp * L
                particles.y[:npart] += yp * L
                particles.s[:npart] += L
                particles.zeta[:npart] += dzeta*L

            track(self, particles, npart)

            # Drift inactive back
            L = self.inactive_back
            if L > 0:
                rpp = particles.rpp[:npart]
                xp = particles.px[:npart] * rpp
                yp = particles.py[:npart] * rpp
                dzeta = particles.rvv[:npart] - ( 1 + ( xp*xp + yp*yp ) / 2 )
                particles.x[:npart] += xp * L
                particles.y[:npart] += yp * L
                particles.s[:npart] += L
                particles.zeta[:npart] += dzeta*L

            particles.reorganize()

