# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

def generate_pencil_on_collimator(line, collimator, num_particles, *, side='+-', impact_parameter=0, 
                                  pencil_spread=1e-6, transverse_impact_parameter=0.,
                                  transverse_spread_sigma=0.01, longitudinal=None,
                                  longitudinal_betatron_cut=None, sigma_z=7.61e-2):
#     if not self.openings_set:
#         raise ValueError("Need to set collimator openings before generating pencil distribution!")
#     if not self.tracker_ready:
#         raise Exception("Please build tracker before generating pencil distribution!")
#     if transverse_impact_parameter != 0.:
#         raise NotImplementedError

    # TODO: check collimator in colldb and installed!

    if collimator.side == 'left':
        side = '+'
    if collimator.side == 'right':
        side = '-'

    # Define the plane
    angle = collimator.angle
    if abs(np.mod(angle-90,180)-90) < 1e-6:
        plane = 'x'
    elif abs(np.mod(angle,180)-90) < 1e-6:
        plane = 'y'
    else:
        raise NotImplementedError("Pencil beam on a skew collimator not yet supported!")

    nemitt_x   = self.colldb.emittance[0]
    nemitt_y   = self.colldb.emittance[1]

    # Is it converging or diverging?
    is_converging = self.colldb._optics.loc[self.s_active_front[collimator], 'alf' + plane ] > 0
    print(f"Collimator {collimator} is {'con' if is_converging else 'di'}verging.")
    if is_converging:
        # pencil at front of jaw
        match_at_s = self.s_active_front[collimator]
        sigma      = self.colldb._beam_size_front[collimator]
    else:
        # pencil at back of jaw
        match_at_s = self.s_active_back[collimator]
        sigma      = self.colldb._beam_size_back[collimator]

    dr_sigmas = pencil_spread/sigma

    # Generate 4D coordinates
    if side == '+-':
        num_plus = int(num_particles/2)
        num_min  = int(num_particles - num_plus)
        coords_plus = self._generate_4D_pencil_one_jaw(num_plus, collimator, plane, '+',
                                                       impact_parameter, dr_sigmas,
                                                       transverse_spread_sigma, match_at_s)
        coords_min  = self._generate_4D_pencil_one_jaw(num_min, collimator, plane, '-',
                                                       impact_parameter, dr_sigmas,
                                                       transverse_spread_sigma, match_at_s)
        coords      = [ [*c_plus, *c_min] for c_plus, c_min in zip(coords_plus, coords_min)]
    else:
        coords      = self._generate_4D_pencil_one_jaw(num_particles, collimator, plane, side,
                                                       impact_parameter, dr_sigmas,
                                                       transverse_spread_sigma, match_at_s)
    pencil            = coords[0]
    p_pencil          = coords[1]
    transverse_norm   = coords[2]
    p_transverse_norm = coords[3]

    # Longitudinal plane
    # TODO: make this more general, make this better
    if longitudinal is None:
        delta = 0
        zeta  = 0
    elif longitudinal == 'matched_dispersion':
        if longitudinal_betatron_cut is None:
            cut = 0
        else:
            cut = np.random.uniform(-longitudinal_betatron_cut, longitudinal_betatron_cut, num_particles)
        delta = self.generate_delta_from_dispersion(at_element=collimator, plane=plane, position_mm=pencil,
                                                    betatron_cut=cut)
        zeta  = 0
    elif longitudinal == 'bucket':
        zeta, delta = xp.generate_longitudinal_coordinates(
                num_particles=num_particles, distribution='gaussian', sigma_z=sigma_z, line=self.line
        )
    elif not hasattr(longitudinal, '__iter__'):
        raise ValueError
    elif len(longitudinal) != 2:
        raise ValueError
    elif isinstance(longitudinal, str):
        raise ValueError
    elif isinstance(longitudinal, dict):
        zeta = longitudinal['zeta']
        delta = longitudinal['delta']
    else:
        zeta = longitudinal[0]
        delta = longitudinal[1]

    # Build the particles
    if plane == 'x':
        part = xp.build_particles(
                x=pencil, px=p_pencil, y_norm=transverse_norm, py_norm=p_transverse_norm,
                zeta=zeta, delta=delta, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                line=self.line, at_element=collimator, match_at_s=match_at_s,
                _context=self._buffer.context
        )
    else:
        part = xp.build_particles(
                x_norm=transverse_norm, px_norm=p_transverse_norm, y=pencil, py=p_pencil, 
                zeta=zeta, delta=delta, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                line=self.line, at_element=collimator, match_at_s=match_at_s,
                _context=self._buffer.context
        )

    part._init_random_number_generator()

    return part


def generate_delta_from_dispersion(self, at_element, *, plane, position_mm, betatron_cut=0):
    line = self.line
    if line.tracker is None:
        raise ValueError("Need to build tracker first!")
    if not hasattr(betatron_cut, '__iter__'):
        if hasattr(position_mm, '__iter__'):
            betatron_cut = np.full_like(position_mm, betatron_cut)
    elif not hasattr(position_mm, '__iter__'):
        position_mm = np.full_like(betatron_cut, position_mm)
    elif len(position_mm) != len(betatron_cut):
        raise ValueError
    tw = line.twiss()
    if isinstance(at_element, str):
        idx = line.element_names.index(at_element)
    betagamma = line.particle_ref.beta0[0] * line.particle_ref.gamma0[0]
    if plane == 'x':
        sigma = np.sqrt(tw.betx[idx]*self.colldb.emittance[0]/betagamma)
        delta = (position_mm - betatron_cut*sigma - tw.x[idx]) / tw.dx[idx]
    elif plane == 'y':
        sigma = np.sqrt(tw.bety[idx]*self.colldb.emittance[1]/betagamma)
        delta = (position_mm - betatron_cut*sigma - tw.y[idx]) / tw.dy[idx]
    else:
        raise ValueError("The variable 'plane' needs to be either 'x' or 'y'!")
    return delta
