import numpy as np
from .k2.materials import Material, CrystalMaterial
from .collimators import BaseCollimator

import xobjects as xo


class K2Engine(xo.HybridClass):

    _xofields = {
        'n_alloc':                xo.Int64,
        'random_generator_seed':  xo.Int64,
#         'random_freeze_state':    xo.Int64,  # to be implemented; number of randoms already sampled, such that this can be taken up again later
#         'collimators':            K2Collimator[:],  # Does not work, need a pointer of sorts
    }

    def __init__(self, **kwargs):
        kwargs.setdefault('_n_alloc', 50000)
        kwargs.setdefault('random_generator_seed', np.random.randint(1, 10000000))
#         kwargs.setdefault('random_freeze_state', -1)
#         kwargs.setdefault('collimators', [])
        super().__init__(**kwargs)
        try:
            import xcoll.beam_elements.pyk2 as pyk2
        except ImportError:
            print("Warning: Failed importing pyK2 (did you compile?). " \
                  + "K2collimators will be installed but are not trackable.")
        else:
            pyk2.pyk2_init(n_alloc=self.n_alloc, random_generator_seed=self.random_generator_seed)



# TODO: remove dx, dy, offset, tilt, as this should only be in colldb (and here only the jaw positions)
class K2Collimator(BaseCollimator):
    _xofields = BaseCollimator._xofields | {
        'dpx':        xo.Float64,
        'dpy':        xo.Float64,
        'offset':     xo.Float64,
        'onesided':   xo.Int8,
        'tilt':       xo.Float64[:],  # TODO: how to limit this to length 2
        'material':   Material,
        'k2engine':   K2Engine
    }

    _skip_in_to_dict       = BaseCollimator._skip_in_to_dict
    _store_in_to_dict      = BaseCollimator._store_in_to_dict
    _internal_record_class = BaseCollimator._internal_record_class

    iscollective = True # TODO: will be set to False when fully in C

    def __init__(self, **kwargs):
        kwargs.setdefault('k2engine', K2Engine())
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
#         self.k2engine.collimators += self
        self._reset_random_seed = False


    def reset_random_seed(self):
        self._reset_random_seed = True


    def track(self, particles):  # TODO: write impacts
        npart = particles._num_active_particles
        if npart > self.k2engine.n_alloc:
            raise ValueError(f"Tracking {npart} particles but only {self.k2engine.n_alloc} allocated!")
        if npart == 0:
            return
        if self._reset_random_seed == True:
            reset_seed = self.k2engine.random_generator_seed
            self._reset_random_seed = False
        else:
            reset_seed = -1
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

            from .k2 import track
            track(self, particles, npart, reset_seed)

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
            
