import numpy as np

import xobjects as xo
import xtrack as xt

from ..tables import CollimatorImpacts

# class MetaCollimator(xt.base_element.MetaBeamElement, ABCMeta):
#     pass

class BaseCollimator(xt.BeamElement):#, metaclass=MetaCollimator):
    _xofields = {
        'inactive_front': xo.Float64,
        'active_length': xo.Float64,
        'inactive_back': xo.Float64,
        'jaw_F_L': xo.Float64,
        'jaw_F_R': xo.Float64,
        'jaw_B_L': xo.Float64,
        'jaw_B_R': xo.Float64,
        'dx': xo.Float64,
        'dy': xo.Float64,
        'cos_z': xo.Float64,
        'sin_z': xo.Float64,
        '_active': xo.Int8
    }

    isthick = True
    behaves_like_drift = True
    
    _skip_in_to_dict  = ['_active', 'cos_z', 'sin_z']
    _store_in_to_dict = ['is_active', 'angle']
    _internal_record_class = CollimatorImpacts

    # BaseCollimator should not be used!
    # We make sure this cannot be accidentally done by killing all particles
    _extra_c_sources = ["""
        void BaseCollimator_track_local_particle(BaseCollimatorData el, LocalParticle* part0) {
            //start_per_particle_block (part0->part)
                LocalParticle_set_x(part, 1e9);
                LocalParticle_set_px(part, 1e9);
                LocalParticle_set_y(part, 1e9);
                LocalParticle_set_py(part, 1e9);
                LocalParticle_set_zeta(part, 1e9);
                LocalParticle_update_delta(part, -1);  // zero energy
            //end_per_particle_block (part0->part)
        }
    """]

    def __init__(self, **kwargs):
        # TODO: quick hack to avoid instantiation; did not manage to get it to work correclty with ABC
        if self.__class__.__name__ == 'BaseCollimator':
            raise Exception("Abstract class 'BaseCollimator' cannot be instantiated!")
        kwargs.setdefault('jaw_F_L', 1)
        kwargs.setdefault('jaw_F_R', -1)
        kwargs.setdefault('jaw_B_L', 1)
        kwargs.setdefault('jaw_B_R', -1)
        kwargs.setdefault('inactive_front', 0)
        kwargs.setdefault('inactive_back', 0)
        kwargs.setdefault('dx', 0)
        kwargs.setdefault('dy', 0)
        angle = kwargs.pop('angle', 0)
        anglerad = angle / 180. * np.pi
        kwargs['cos_z'] = np.cos(anglerad)
        kwargs['sin_z'] = np.sin(anglerad)
        is_active = kwargs.pop('is_active', True)
        is_active = 1 if is_active == True else is_active
        is_active = 0 if is_active == False else is_active
        kwargs['_active'] = is_active
        super().__init__(**kwargs)


    @property
    def angle(self):
        return np.arctan2(self.sin_z, self.cos_z) * (180.0 / np.pi)

    @angle.setter
    def angle(self, angle):
        anglerad = angle / 180. * np.pi
        self.cos_z = np.cos(anglerad)
        self.sin_z = np.sin(anglerad)

    @property
    def is_active(self):
        return True if self._active == 1 else False

    @is_active.setter
    def is_active(self, is_active):
        is_active = 1 if is_active == True else is_active
        is_active = 0 if is_active == False else is_active
        self._active = is_active
        if is_active <= 0:
            self.jaw_F_L = 1
            self.jaw_F_R = -1
            self.jaw_B_L = 1
            self.jaw_B_R = -1 

    @property
    def length(self):
        return (self.inactive_front + self.active_length + self.inactive_back)

    # TODO: In principle we are not allowed to backtrack through a collimator
    #       However, the loss refinement will fail if this function is not provided
    def get_backtrack_element(self, _context=None, _buffer=None, _offset=None):
        return xt.Drift(length=-self.length, _context=_context, _buffer=_buffer, _offset=_offset)


