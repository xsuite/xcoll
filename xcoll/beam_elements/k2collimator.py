import numpy as np


class K2Engine:

    def __init__(self, n_alloc, random_generator_seed=None):

        self._n_alloc = n_alloc

        if random_generator_seed is None:
            self.random_generator_seed = np.random.randint(1, 100000)
        else:
            self.random_generator_seed = random_generator_seed

        try:
            import xcoll.beam_elements.pyk2 as pyk2
        except ImportError:
            print("Warning: Failed importing pyK2 (did you compile?)." \
                  + "K2collimators will be installed but are not trackable.")
        else:
            pyk2.pyk2_init(n_alloc=n_alloc, random_generator_seed=self.random_generator_seed)

    @property
    def n_alloc(self):
        return self._n_alloc


class K2Collimator:

    iscollective = True
    isthick = True
    behaves_like_drift = True

    def __init__(self, *, k2engine, active_length, angle, inactive_front=0, inactive_back=0,
                 jaw_F_L=1, jaw_F_R=-1, jaw_B_L=1, jaw_B_R=-1, onesided=False, dx=0, dy=0, dpx=0, dpy=0, offset=0, tilt=[0,0],
                 is_active=True, record_impacts=False, material=None):

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
        self.record_impacts = record_impacts
        self._reset_random_seed = False

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
        import xcoll.beam_elements.pyk2 as pyk2
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


            # =================================================================== #
            # ===============================  K2  ============================== #
            # =================================================================== #
            #  Even though K2 has info on the length of the collimator, and even  #
            #  uses that info, it is still a thin element in the original code.   #
            #  Hence, our thick tracking will always be equivalent to:            #
            #    - drift first half length                                        #
            #    - scatter                                                        #
            #    - drift second half length                                       #
            # =================================================================== #

            # Go to collimator reference system (subtract closed orbit)
            # For standard K2, the CO and gap should be taken in the centre
            x_part  = particles.x[:npart] - self.dx
            xp_part = (particles.px[:npart] - self.dpx) * particles.rpp[:npart]
            y_part  = particles.y[:npart] - self.dy
            yp_part = (particles.py[:npart] - self.dpy) * particles.rpp[:npart]
            s_part  = 0 * x_part
            e_part  = particles.energy[:npart].copy() / 1e9 # Energy in GeV
            rpp_in  = particles.rpp[:npart].copy()
            rvv_in  = particles.rvv[:npart].copy()

            # Initialise arrays for FORTRAN call
            part_hit = np.zeros(len(x_part), dtype=np.int32)
            part_abs = np.zeros(len(x_part), dtype=np.int32)
            part_impact = np.zeros(len(x_part), dtype=np.float)
            part_indiv = np.zeros(len(x_part), dtype=np.float)
            part_linteract = np.zeros(len(x_part), dtype=np.float)
            nhit_stage = np.zeros(len(x_part), dtype=np.int32)
            nabs_type = np.zeros(len(x_part), dtype=np.int32)
            linside = np.zeros(len(x_part), dtype=np.int32)

            if self.jaw_F_L != self.jaw_B_L or self.jaw_F_R != self.jaw_B_R:
                raise NotImplementedError
            opening = self.jaw_F_L - self.jaw_F_R
            offset = self.offset + ( self.jaw_F_L + self.jaw_F_R )/2

            matID = pyk2.materials[self.material]['ID']
            
            pyk2.pyk2_run(x_particles=x_part,
                      xp_particles=xp_part,
                      y_particles=y_part,
                      yp_particles=yp_part,
                      s_particles=s_part,
                      p_particles=e_part,              # confusing: this is ENERGY not momentum
                      part_hit=part_hit,
                      part_abs=part_abs,
                      part_impact=part_impact,         # impact parameter
                      part_indiv=part_indiv,           # particle divergence
                      part_linteract=part_linteract,   # interaction length
                      nhit_stage=nhit_stage,
                      nabs_type=nabs_type,
                      linside=linside,
                      matid=matID,
                      is_crystal=False,
                      c_length=self.active_length,
                      c_rotation=self.angle/180.*np.pi,
                      c_aperture=opening,
                      c_offset=offset,
                      c_tilt=self.tilt,
                      c_enom=particles.energy0[0]/1e6, # Reference energy
                      onesided=self.onesided,
                      random_generator_seed=reset_seed,
                      )

            # Masks of hit and survived particles
            lost = part_abs > 0
            hit = part_hit > 0
            not_hit = ~hit
            not_lost = ~lost
            survived_hit = hit & (~lost)
            
#             print(lost[:40])
#             print(hit[:40])
#             print(survived_hit[:40])
#             print(part_impact[:40])
#             print(part_indiv[:40])
#             print(part_linteract[:40])  #  This is how much of the collimator it traversed -- correct? Or only what was left?
#             print(nabs_type[:40])

            # Update energy    ---------------------------------------------------
            # Only particles that hit the jaw and survived need to be updated
            ptau_out = particles.ptau[:npart]
            e0 = particles.energy0[:npart]
            beta0 = particles.beta0[:npart]
            ptau_out[survived_hit] = (
                    e_part[survived_hit] * 1e9 - e0[survived_hit]
                ) / (
                    e0[survived_hit] * beta0[survived_hit]
                )
            particles.ptau[:npart] = ptau_out

#             # Sanity check: non-hit particles get no angle correction:
#             print(np.allclose(rvv_in[not_hit], particles.rvv[not_hit], atol=1e-15, rtol=0))
#             print(np.allclose(rpp_in[not_hit], particles.rpp[not_hit], atol=1e-15, rtol=0))

            # Update momenta    --------------------------------------------------
            # Absorbed particles get coordinates set to the entrance of collimator
            px_out = particles.px[:npart].copy()
            py_out = particles.py[:npart].copy()
            # Non-hit particles are just drifting
            px_out[not_hit] = xp_part[not_hit]/particles.rpp[not_hit]
            py_out[not_hit] = yp_part[not_hit]/particles.rpp[not_hit]
            # Hit and survived particles need correcting:
            # Note that K2 did not update angles yet with new energy!
            # So need to do xp' = xp * p_in / p_out
            # (see collimation.f90 line 1709 and mod_particles.f90 line 210)
            # and then transform them to px = xp' * p_out/p0 = xp * p_in/p0 = xp / rpp_in
            px_out[survived_hit] = xp_part[survived_hit]/rpp_in[survived_hit]
            py_out[survived_hit] = yp_part[survived_hit]/rpp_in[survived_hit]

            # Update and correct transversal coordinates    ----------------------
            # Absorbed particles get coordinates set to the entrance of collimator
            x_out = particles.x[:npart].copy()
            y_out = particles.y[:npart].copy()
            # Non-hit particles are just drifting
            x_out[not_hit] = x_part[not_hit]
            y_out[not_hit] = y_part[not_hit]
            # Hit and survived particles need correcting:
            # In SixTrack K2, particles are first backtracked to the centre of the collimator,
            # then the angles are updated, then they are shifted back to the end of the collimator
            # 1) Backtrack to centre with non-updated angle: x += - xp L/2
            # 2) Track to end with updated angle: x += xp' L/2
            # Total: x += xp * L/2 * ( p_in / p_out - 1 )
            #        x += xp * L/2 * ( rpp / rpp_in - 1 )
            correction = self.active_length/2*(particles.rpp/rpp_in-1)
            x_out[survived_hit] = x_part[survived_hit] + xp_part[survived_hit]*correction[survived_hit]
            y_out[survived_hit] = y_part[survived_hit] + yp_part[survived_hit]*correction[survived_hit]
            
            # Move back from collimator frame
            x_out += self.dx
            y_out += self.dy
            px_out += self.dpx
            py_out += self.dpy
            particles.x[:npart] = x_out
            particles.y[:npart] = y_out
            particles.px[:npart] = px_out
            particles.py[:npart] = py_out

            # Update longitudinal coordinate zeta    -----------------------------
            # Absorbed particles get coordinates set to the entrance of collimator
            rvv_out = particles.rvv[:npart]
            zeta_out = particles.zeta[:npart]
            # Non-hit particles are just drifting (not yet drifted in K2, so do here)
            zeta_out[not_hit] += self.active_length*(
                                rvv_out[not_hit] - (1 + ( xp_part[not_hit]**2 + yp_part[not_hit]**2)/2 ) 
                            )
            # Hit and survived particles need correcting:
            # First we drift half the length with the old angles, then half the length
            # with the new angles
            zeta_out[survived_hit] += self.active_length/2*( rvv_in[survived_hit] - 
                                        (1 + (xp_part[survived_hit]**2 + yp_part[survived_hit]**2)/2) 
                                      )
            correction = particles.rpp[survived_hit]/rpp_in[survived_hit]
            zeta_out[survived_hit] += self.active_length/2*( rvv_out[survived_hit] -
                                        (1 + (xp_part[survived_hit]**2 + yp_part[survived_hit]**2)/2*correction**2)
                                      )
            # update
            particles.zeta[:npart] = zeta_out

            # Update s    --------------------------------------------------------
            s_out = particles.s[:npart]
            s_out[not_lost] += self.length
            particles.s[:npart] = s_out

            # Update state    ----------------------------------------------------
            state_out = particles.state[:npart].copy()
            state_out[lost] = -333
            particles.state[:npart] = state_out

            # TODO update records

            # =================================================================== #


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
        thisdict = {}
        thisdict['__class__'] = 'K2Collimator'
        thisdict['n_alloc'] = self._k2engine.n_alloc
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
        thisdict['record_impacts'] = 1 if self.record_impacts else 0
        thisdict['material'] = self.material
        return thisdict


    @classmethod
    def from_dict(cls, thisdict, *, engine=None):
        if engine is None:
            print("Warning: no engine given! Creating a new one...")
            engine = K2Engine(thisdict['n_alloc'], thisdict['random_generator_seed'])
        else:
            if engine.n_alloc != thisdict['n_alloc']:
                raise ValueError("The provided engine is incompatible with the engine of the "\
                                 "stored element: n_alloc is different.")
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
            record_impacts = True if thisdict['record_impacts']==1 else False,
            material = thisdict['material']
        )
            
