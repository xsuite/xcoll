import numpy as np
from .k2 import track
from .k2 import materials

class K2Engine:

    def __init__(self, random_generator_seed=None):

        if random_generator_seed is None:
            self.random_generator_seed = np.random.randint(1, 100000)
        else:
            self.random_generator_seed = random_generator_seed

        try:
            import xcoll.beam_elements.pyk2 as pyk2
        except ImportError:
            print("Warning: Failed importing pyK2 (did you compile?). " \
                  + "K2collimators will be installed but are not trackable.")
        else:
            pyk2.pyk2_init(random_generator_seed=self.random_generator_seed)
        import xcoll.beam_elements.k2.k2_random as kr
        kr.__index__ = -1


class K2Collimator:

    iscollective = True
    isthick = True
    behaves_like_drift = True

    def __init__(self, *, k2engine, active_length, angle, inactive_front=0, inactive_back=0,
                 jaw_F_L=1, jaw_F_R=-1, jaw_B_L=1, jaw_B_R=-1, onesided=False, dx=0, dy=0, dpx=0, dpy=0, offset=0, tilt=[0,0],
                 is_active=True, impacts=None, material=None):

        self._k2engine = k2engine
        self.active_length = active_length
        self.inactive_front = inactive_front
        self.inactive_back = inactive_back
        self.angle = angle
        self.jaw_F_L = jaw_F_L
        self.jaw_F_R = jaw_F_R
        self.jaw_B_L = jaw_B_L
        self.jaw_B_R = jaw_B_R
        self.onesided = onesided
        self.dx = dx
        self.dy = dy
        self.dpx = dpx
        self.dpy = dpy
        self.offset = offset
        if hasattr(tilt, '__iter__'):
            if isinstance(tilt, str):
                raise ValueError("Variable tilt has to be a number or array of numbers!")
            elif len(tilt) == 1:
                tilt = [tilt[0], tilt[0]]
        else:
            tilt = [tilt, tilt]
        self.tilt = np.array(tilt, dtype=np.float64)
        self.is_active = is_active
        self.material = material
        self.impacts = impacts
        self._reset_random_seed = False

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if not material in materials.keys():
            raise ValueError(f"Argument '{material}' is not in the materials database!")
        self._material = material

    @property
    def impacts(self):
        return self._impacts

    @property
    def record_impacts(self):
        return self._record_impacts

    @impacts.setter
    def impacts(self, impacts):
        self._record_impacts = False if impacts is None else True
        self._impacts = impacts

    @property
    def k2engine(self):
        return self._k2engine

    @property
    def is_active(self):
        return self._active

    @is_active.setter
    def is_active(self, is_active):
        self._active = is_active
        if not is_active:
            self.jaw_F_L = 1
            self.jaw_F_R = -1
            self.jaw_B_L = 1
            self.jaw_B_R = -1

    def reset_random_seed(self):
        self._reset_random_seed = True

    @property
    def length(self):
        return self.active_length + self.inactive_front + self.inactive_back

    def track(self, particles):  # TODO: write impacts
        npart = particles._num_active_particles
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
        
    def to_dict(self):
        # TODO how to save ref to impacts?
        thisdict = {}
        thisdict['__class__'] = 'K2Collimator'
        thisdict['random_generator_seed'] = self._k2engine.random_generator_seed
        thisdict['active_length'] = self.active_length
        thisdict['inactive_front'] = self.inactive_front
        thisdict['inactive_back'] = self.inactive_back
        thisdict['angle'] = self.angle
        thisdict['jaw_F_L'] = self.jaw_F_L
        thisdict['jaw_F_R'] = self.jaw_F_R
        thisdict['jaw_B_L'] = self.jaw_B_L
        thisdict['jaw_B_R'] = self.jaw_B_R
        thisdict['onesided'] = 1 if self.onesided else 0
        thisdict['dx'] = self.dx
        thisdict['dy'] = self.dy
        thisdict['dpx'] = self.dpx
        thisdict['dpy'] = self.dpy
        thisdict['offset'] = self.offset
        thisdict['tilt'] = self.tilt
        thisdict['is_active'] = 1 if self.is_active else 0
        thisdict['material'] = self.material
        return thisdict


    @classmethod
    def from_dict(cls, thisdict, *, engine=None):
         # TODO how to get ref to impacts?
        if engine is None:
            print("Warning: no engine given! Creating a new one...")
            engine = K2Engine(thisdict['random_generator_seed'])
        else:
            if engine.random_generator_seed != thisdict['random_generator_seed']:
                raise ValueError("The provided engine is incompatible with the engine of the "\
                                 "stored element: random_generator_seed is different.")
        return cls(
            k2engine = engine,
            active_length = thisdict['active_length'],
            inactive_front = thisdict['inactive_front'],
            inactive_back = thisdict['inactive_back'],
            angle = thisdict['angle'],
            jaw_F_L = thisdict['jaw_F_L'],
            jaw_F_R = thisdict['jaw_F_R'],
            jaw_B_L = thisdict['jaw_B_L'],
            jaw_B_R = thisdict['jaw_B_R'],
            onesided = True if thisdict['onesided']==1 else False,
            dx = thisdict['dx'],
            dy = thisdict['dy'],
            dpx = thisdict['dpx'],
            dpy = thisdict['dpy'],
            offset = thisdict['offset'],
            tilt = thisdict['tilt'],
            is_active = True if thisdict['is_active']==1 else False,
            material = thisdict['material']
        )


class K2Crystal:

    iscollective = True
    isthick = True
    behaves_like_drift = True

    def __init__(self, *, k2engine, active_length, angle, align_angle, inactive_front=0, inactive_back=0,
                 jaw_F_L=1, jaw_F_R=-1, jaw_B_L=1, jaw_B_R=-1, onesided=False, dx=0, dy=0, dpx=0, dpy=0, offset=0, tilt=[0,0],
                 bend=0, xdim=0, ydim=0, thick=0, crytilt=0, miscut=0, orient=0,
                 is_active=True, impacts=None, material=None):

        self._k2engine = k2engine
        self.active_length = active_length
        self.inactive_front = inactive_front
        self.inactive_back = inactive_back
        self.angle = angle
        self.jaw_F_L = jaw_F_L
        self.jaw_F_R = jaw_F_R
        self.jaw_B_L = jaw_B_L
        self.jaw_B_R = jaw_B_R
        self.onesided = onesided
        self.dx = dx
        self.dy = dy
        self.dpx = dpx
        self.dpy = dpy
        self.offset = offset
        if hasattr(tilt, '__iter__'):
            if isinstance(tilt, str):
                raise ValueError("Variable tilt has to be a number or array of numbers!")
            elif len(tilt) == 1:
                tilt = [tilt[0], tilt[0]]
        else:
            tilt = [tilt, tilt]
        self.tilt = np.array(tilt, dtype=np.float64)
        self.is_active = is_active
        self.material = material
        self.bend = bend
        self.xdim = xdim
        self.ydim = ydim
        self.thick = thick
        self.crytilt = crytilt
        self.align_angle = align_angle    #  = - sqrt(eps/beta)*alpha*nsigma
        self.miscut = miscut
        self.orient = orient
        self.impacts = impacts
        self._reset_random_seed = False

    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material):
        if not material in materials.keys():
            raise ValueError(f"Argument '{material}' is not in the materials database!")
        if not materials[material]['can_be_crystal']:
            raise ValueError(f"Material '{material}' cannot be used for crystals!")
        self._material = material
    
    @property
    def impacts(self):
        return self._impacts

    @property
    def record_impacts(self):
        return self._record_impacts

    @impacts.setter
    def impacts(self, impacts):
        self._record_impacts = False if impacts is None else True
        self._impacts = impacts

    @property
    def k2engine(self):
        return self._k2engine

    @property
    def is_active(self):
        return self._active

    @is_active.setter
    def is_active(self, is_active):
        self._active = is_active
        if not is_active:
            self.jaw_F_L = 1
            self.jaw_F_R = -1
            self.jaw_B_L = 1
            self.jaw_B_R = -1

    def reset_random_seed(self):
        self._reset_random_seed = True

    @property
    def length(self):
        return self.active_length + self.inactive_front + self.inactive_back

    def track(self, particles):  # TODO: write impacts
        npart = particles._num_active_particles
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
        
    def to_dict(self):
        # TODO how to save ref to impacts?
        thisdict = {}
        thisdict['__class__'] = 'K2Crystal'
        thisdict['random_generator_seed'] = self._k2engine.random_generator_seed
        thisdict['active_length'] = self.active_length
        thisdict['inactive_front'] = self.inactive_front
        thisdict['inactive_back'] = self.inactive_back
        thisdict['angle'] = self.angle
        thisdict['jaw_F_L'] = self.jaw_F_L
        thisdict['jaw_F_R'] = self.jaw_F_R
        thisdict['jaw_B_L'] = self.jaw_B_L
        thisdict['jaw_B_R'] = self.jaw_B_R
        thisdict['onesided'] = 1 if self.onesided else 0
        thisdict['dx'] = self.dx
        thisdict['dy'] = self.dy
        thisdict['dpx'] = self.dpx
        thisdict['dpy'] = self.dpy
        thisdict['offset'] = self.offset
        thisdict['tilt'] = self.tilt
        thisdict['is_active'] = 1 if self.is_active else 0
        thisdict['material'] = self.material
        thisdict['bend'] = self.bend
        thisdict['xdim'] = self.xdim
        thisdict['ydim'] = self.ydim
        thisdict['thick'] = self.thick
        thisdict['crytilt'] = self.crytilt
        thisdict['miscut'] = self.miscut
        thisdict['orient'] = self.orient
        thisdict['align_angle'] = self.align_angle

        return thisdict


    @classmethod
    def from_dict(cls, thisdict, *, engine=None):
         # TODO how to get ref to impacts?
        if engine is None:
            print("Warning: no engine given! Creating a new one...")
            engine = K2Engine(thisdict['random_generator_seed'])
        else:
            if engine.random_generator_seed != thisdict['random_generator_seed']:
                raise ValueError("The provided engine is incompatible with the engine of the "\
                                 "stored element: random_generator_seed is different.")

        return cls(
            k2engine = engine,
            active_length = thisdict['active_length'],
            inactive_front = thisdict['inactive_front'],
            inactive_back = thisdict['inactive_back'],
            angle = thisdict['angle'],
            jaw_F_L = thisdict['jaw_F_L'],
            jaw_F_R = thisdict['jaw_F_R'],
            jaw_B_L = thisdict['jaw_B_L'],
            jaw_B_R = thisdict['jaw_B_R'],
            onesided = True if thisdict['onesided']==1 else False,
            dx = thisdict['dx'],
            dy = thisdict['dy'],
            dpx = thisdict['dpx'],
            dpy = thisdict['dpy'],
            offset = thisdict['offset'],
            tilt = thisdict['tilt'],
            is_active = True if thisdict['is_active']==1 else False,
            material = thisdict['material'],
            bend = thisdict['bend'],
            xdim = thisdict['xdim'],
            ydim = thisdict['ydim'],
            thick = thisdict['thick'],
            crytilt = thisdict['crytilt'],
            miscut = thisdict['miscut'],
            orient = thisdict['orient'],
            align_angle = thisdict['align_angle']
        )
