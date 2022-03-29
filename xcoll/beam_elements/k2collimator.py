import numpy as np

class K2Engine:

    def __init__(self, n_alloc, colldb_filename, random_generator_seed=None):

        self._n_alloc = n_alloc
        self._colldb_filename = colldb_filename

        if random_generator_seed is None:
            self.random_generator_seed = np.random.randint(1, 100000)
        else:
            self.random_generator_seed = random_generator_seed

        import pyk2
        pyk2.pyk2_init(n_alloc=n_alloc, colldb_input_fname=colldb_filename,
               random_generator_seed=self.random_generator_seed)
        pyk2._active_engine = self
        
        
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

    @property
    def colldb_filename(self):
        return self._colldb_filename


class K2Collimator:

    iscollective = True
    isthick = True
#     behaves_like_drift = True

    def __init__(self, *, k2engine, icoll, active_length, inactive_front, inactive_back, angle,
                 jaw, onesided=False, dx=0, dy=0, dpx=0, dpy=0, offset=0, tilt=None):

        self._k2engine = k2engine
        self.icoll = icoll
        self.active_length = active_length
        self.inactive_front = inactive_front
        self.inactive_back = inactive_back
        self.angle = angle
        self.jaw = jaw
        self.onesided = onesided
        self.dx = dx
        self.dy = dy
        self.dpx = dpx
        self.dpy = dpy
        self.offset = offset
        self.tilt = tilt

    @property
    def k2engine(self):
        return self._k2engine

    @property
    def length(self):
        return self.active_length + self.inactive_front + self.inactive_back

    def track(self, particles):

        import pyk2
        if pyk2._active_engine is not self.k2engine:
            raise ValueError(f"Collimator has different K2Engine than the main initiated one!")
        npart = particles._num_active_particles
        if npart > self.k2engine.n_alloc:
            raise ValueError(f"Tracking {npart} particles but only {self.k2engine.n_alloc} allocated!")

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

        # Go to collimator reference system (subtract closed orbit)
        x_part = particles.x[:npart] - self.dx
        xp_part = (particles.px[:npart] - self.dpx) * particles.rpp[:npart]
        y_part = particles.y[:npart] - self.dy
        yp_part = (particles.py[:npart] - self.dpy) * particles.rpp[:npart]
        s_part = 0 * x_part
        p_part = particles.energy[:npart] / 1e9 # Energy (not momentum) in GeV

        # Initialise arrays for FORTRAN call
        part_hit_pos = np.zeros(len(x_part), dtype=np.int32)
        part_hit_turn = np.zeros(len(x_part), dtype=np.int32)
        part_abs_pos = np.zeros(len(x_part), dtype=np.int32)
        part_abs_turn = np.zeros(len(x_part), dtype=np.int32)
        part_impact = np.zeros(len(x_part), dtype=np.float)
        part_indiv = np.zeros(len(x_part), dtype=np.float)
        part_linteract = np.zeros(len(x_part), dtype=np.float)
        nhit_stage = np.zeros(len(x_part), dtype=np.int32)
        nabs_type = np.zeros(len(x_part), dtype=np.int32)
        linside = np.zeros(len(x_part), dtype=np.int32)

        # `linside` is an array of logicals in fortran. Beware of the fortran converion:
        # True <=> -1 (https://stackoverflow.com/questions/39454349/numerical-equivalent-of-true-is-1)

        pyk2.pyk2_run(x_particles=x_part,
                      xp_particles=xp_part,
                      y_particles=y_part,
                      yp_particles=yp_part,
                      s_particles=s_part,
                      p_particles=p_part,              # confusing: this is ENERGY not momentum
                      part_hit_pos=part_hit_pos,       # ignore: sixtrack element of impact
                      part_hit_turn=part_hit_turn,     # ignore: turn of impact
                      part_abs_pos=part_abs_pos,       # ignore: sixtrack element of absorption
                      part_abs_turn=part_abs_turn,     # ignore: turn of absorption
                      part_impact=part_impact,         # impact parameter
                      part_indiv=part_indiv,           # particle divergence
                      part_linteract=part_linteract,   # interaction length
                      nhit_stage=nhit_stage,
                      nabs_type=nabs_type,
                      linside=linside,
                      icoll=self.icoll,
                      ie=1,                            # ignore: structure element index
                      c_length=self.active_length,
                      c_rotation=self.angle/180.*np.pi,
                      c_aperture=2*self.jaw,
                      c_offset=self.offset,
                      c_tilt=np.array([0,0], dtype=np.float64),
                      c_enom=particles.energy0[0]/1e9, # Reference energy
                      onesided=self.onesided,
                      random_generator_seed=-1, # skips rng re-initlization
                      )

        # Masks of hit and survived particles
        mask_lost = part_abs_turn > 0
        mask_hit = part_hit_turn > 0
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
        
        
