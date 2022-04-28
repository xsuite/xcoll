import numpy as np
from pyk2 import k2_track

class K2Engine:

    def __init__(self, n_alloc, random_generator_seed=None):

        self._n_alloc = n_alloc

        if random_generator_seed is None:
            self.random_generator_seed = np.random.randint(1, 100000)
        else:
            self.random_generator_seed = random_generator_seed

        import xcoll.beam_elements.pyk2 as pyk2
        pyk2.pyk2_init(n_alloc=n_alloc, random_generator_seed=self.random_generator_seed)
        
        
# call collmat_init                   ! Set default values for collimator materials
# call cdb_readCollDB                 ! Read the collimator DB
# call cdb_setLHCOnesided(do_oneside) ! Set LHC onesided collimators
# call k2coll_init
# if(coll_hasCrystal) then
#     call cry_init
# end if

    @property
    def n_alloc(self):
        return self._n_alloc


class K2Collimator:

    iscollective = True
    isthick = True
#     behaves_like_drift = True

    def __init__(self, *, k2engine, active_length, inactive_front, inactive_back, angle, is_active=True,
                 jaw_F_L=1, jaw_F_R=-1, jaw_B_L=1, jaw_B_R=-1, onesided=False, dx=0, dy=0, dpx=0, dpy=0, offset=0, tilt=0, material=None):

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
        self.tilt = tilt
        self._active = is_active
        self.material = material
        self.is_crystal = False

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

    @property
    def length(self):
        return self.active_length + self.inactive_front + self.inactive_back

    def track(self, particles):  # TODO: write impacts
        import xcoll.beam_elements.pyk2 as pyk2
        npart = particles._num_active_particles
        if npart > self.k2engine.n_alloc:
            raise ValueError(f"Tracking {npart} particles but only {self.k2engine.n_alloc} allocated!")
        
        if not self._active:
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

            x_part, xp_part, y_part, yp_part, p_part, s_part, part_hit, part_abs = k2_track(
                            material=self.material,
                            particles=particles,
                            closed_orbit=[self.dx,self.dy,self.dpx,self.dpy],
                            jaws=[self.jaw_F_L,self.jaw_F_R,self.jaw_B_L,self.jaw_B_R]
                            offset=self.offset, npart=npart, length=self.active_length,
                            is_crystal=self.is_crystal,
                            onesided=self.onesided
                        )
            
            # Masks of hit and survived particles
            mask_lost = abs > 0
            mask_hit = part_hit > 0
            mask_not_hit = ~mask_hit
            mask_survived_hit = mask_hit & (~mask_lost)

            state_out = particles.state[:npart].copy()
            state_out[mask_lost] = -333
            particles.state[:npart] = state_out

            # Update particle energy
            ptau_out = particles.ptau[:npart].copy()
            e0 = particles.energy0[:npart]
            beta0 = particles.beta0[:npart]
            ptau_out[mask_survived_hit] = (
                    p_part[mask_survived_hit] * 1e9 - e0[mask_survived_hit]
                ) / (
                    e0[mask_survived_hit] * beta0[mask_survived_hit]
                )
            particles.ptau[:npart] = ptau_out

            # Update transversal coordinates (moving back from collimator frame)
            x_part_out = particles.x[:npart].copy()
            x_part_out[~mask_lost] = x_part[~mask_lost] + self.dx
            particles.x[:npart] = x_part_out

            y_part_out = particles.y[:npart].copy()
            y_part_out[~mask_lost] = y_part[~mask_lost] + self.dy
            particles.y[:npart] = y_part_out

            rpp_out = particles.rpp[:npart]
            px_out = particles.px[:npart]
            px_out[mask_survived_hit] = xp_part[mask_survived_hit]/rpp_out[mask_survived_hit]
            particles.px[:npart] = px_out

            py_out = particles.py[:npart]
            py_out[mask_survived_hit] = yp_part[mask_survived_hit]/rpp_out[mask_survived_hit]
            particles.py[:npart] = py_out

            # Update longitudinal coordinate zeta
            # Note: they have NOT been drifted yet in K2, so we have to do that manually
            rvv_out = particles.rvv[:npart]
            zeta_out = particles.zeta[:npart]
            zeta_out[mask_not_hit] += self.length*(
                                rvv_out[mask_not_hit] - (1 + ( xp_part[mask_not_hit]**2 + yp_part[mask_not_hit]**2)/2 ) 
                            )
            particles.zeta[:npart] = zeta_out

            s_out = particles.s[:npart]
            s_out[mask_not_hit] += self.length
            particles.s[:npart] = s_out

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
        
        
