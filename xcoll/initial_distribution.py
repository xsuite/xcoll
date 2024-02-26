# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #
import numpy as np

import xtrack as xt
import xobjects as xo
import xpart as xp

from .beam_elements import _all_collimator_types

def generate_pencil_on_collimator(line, collimator_name, num_particles, *, side='+-', impact_parameter=0, 
                                  pencil_spread=1e-6, transverse_impact_parameter=0.,
                                  transverse_spread_sigma=0.01, longitudinal=None,
                                  longitudinal_betatron_cut=None, sigma_z=7.61e-2,emittance, 
                                  beta_gamma_rel):
    """
    Generate a pencil beam on a collimator. 
    """
    # TODO: horrible solution for collimation_manager.openings_set, fix this later
    
    if line.tracker is None:
        raise Exception("Please build tracker before generating pencil distribution!")
    if transverse_impact_parameter != 0.:
        raise NotImplementedError
    if not isinstance(collimator_name, str):
        raise ValueError("Need to provide a string!")
    else:
        collimator = line[collimator_name]
    if not isinstance(collimator, tuple(_all_collimator_types)):
        raise ValueError("Need to provide a valid collimator!")
    if collimator.active == 0:
        raise Exception("Need to set collimator openings before generating pencil distribution!")
    
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
    
    nemitt_x = emittance[0] 
    nemitt_y = emittance[1]

    # Is it converging or diverging?
    tw = line.twiss()
    s_active_front = line.get_s_position(collimator_name)
    s_active_back  = line.get_s_position(collimator_name) + collimator.active_length
    is_converging  = tw['alfx',collimator_name] > 0
    print(f"Collimator {collimator} is {'con' if is_converging else 'di'}verging.")

    if is_converging:
        # pencil at front of jaw
        print('Im in the front')
        match_at_s = s_active_front
        betax      = tw['betx',collimator_name + '_upstream_aper']
        betay      = tw['bety',collimator_name + '_upstream_aper']
        sigma_x    = np.sqrt(betax*nemitt_x/beta_gamma_rel)
        sigma_y    = np.sqrt(betay*nemitt_y/beta_gamma_rel)
        sigma      = np.sqrt( (sigma_x*np.cos(np.float_(collimator.angle)*np.pi/180))**2
                            + (sigma_y*np.sin(np.float_(collimator.angle)*np.pi/180))**2)
    else:
        # pencil at back of jaw
        print('Im in the back')
        match_at_s = s_active_back
        betax      = tw['betx',collimator_name + '_downstream_aper']
        betay      = tw['bety',collimator_name + '_downstream_aper']
        sigma_x    = np.sqrt(betax*nemitt_x/beta_gamma_rel)
        sigma_y    = np.sqrt(betay*nemitt_y/beta_gamma_rel)
        sigma      = np.sqrt( (sigma_x*np.cos(np.float_(collimator.angle)*np.pi/180))**2
                            + (sigma_y*np.sin(np.float_(collimator.angle)*np.pi/180))**2)

    dr_sigmas = pencil_spread/sigma

    # Generate 4D coordinates
    if side == '+-':
        num_plus = int(num_particles/2)
        num_min  = int(num_particles - num_plus)
        coords_plus = generate_4D_pencil_one_jaw(num_plus, line, collimator_name, plane, '+',
                                                impact_parameter, dr_sigmas,
                                                transverse_spread_sigma, match_at_s,emittance)
        coords_min  = generate_4D_pencil_one_jaw(num_min, line, collimator_name, plane, '-',
                                                impact_parameter, dr_sigmas,
                                                transverse_spread_sigma, match_at_s,emittance)
        coords      = [ [*c_plus, *c_min] for c_plus, c_min in zip(coords_plus, coords_min)]
    else:
        coords      = generate_4D_pencil_one_jaw(num_particles, line, collimator_name, plane, side,
                                                impact_parameter, dr_sigmas,
                                                transverse_spread_sigma, match_at_s,emittance)
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
        delta = generate_delta_from_dispersion(emittance=emittance,at_element=collimator_name, line=line, plane=plane, position_mm=pencil,
                                                betatron_cut=cut)
        zeta  = 0
    elif longitudinal == 'bucket':
        zeta, delta = xp.generate_longitudinal_coordinates(
                num_particles=num_particles, distribution='gaussian', sigma_z=sigma_z, line=line
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
                line=line, at_element=collimator_name, match_at_s=match_at_s,
                _context=collimator._buffer.context
        )
    else:
        part = xp.build_particles(
                x_norm=transverse_norm, px_norm=p_transverse_norm, y=pencil, py=p_pencil, 
                zeta=zeta, delta=delta, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                line=line, at_element=collimator_name, match_at_s=match_at_s,
                _context=collimator._buffer.context
        )

    part._init_random_number_generator()

    return part


def generate_delta_from_dispersion(at_element,line, *, plane, position_mm,emittance, betatron_cut=0):
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
        sigma = np.sqrt(tw.betx[idx]*emittance[0]/betagamma)
        delta = (position_mm - betatron_cut*sigma - tw.x[idx]) / tw.dx[idx]
    elif plane == 'y':
        sigma = np.sqrt(tw.bety[idx]*emittance[1]/betagamma)
        delta = (position_mm - betatron_cut*sigma - tw.y[idx]) / tw.dy[idx]
    else:
        raise ValueError("The variable 'plane' needs to be either 'x' or 'y'!")
    return delta

def generate_4D_pencil_one_jaw(num_particles,line,collimator_name, plane, side, impact_parameter,
                                    dr_sigmas, transverse_spread_sigma, match_at_s,emittance):
        if plane == 'x':
            co_pencil     = line[collimator_name].ref_x
            co_transverse = line[collimator_name].ref_y
        else:
            co_pencil     = line[collimator_name].ref_y
            co_transverse = line[collimator_name].ref_x

        nemitt_x   = emittance[0]
        nemitt_y   = emittance[1]

        if side == '+':
            absolute_cut = line[collimator_name].jaw_LU + co_pencil + impact_parameter
        elif side == '-':
            absolute_cut = line[collimator_name].jaw_RU + co_pencil - impact_parameter

        # Collimator plane: generate pencil distribution
        pencil, p_pencil = xp.generate_2D_pencil_with_absolute_cut(num_particles,
                        plane=plane, absolute_cut=absolute_cut, dr_sigmas=dr_sigmas,
                        side=side, line=line, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                        at_element=collimator_name, match_at_s=match_at_s
        )

        # Other plane: generate gaussian distribution in normalized coordinates
        transverse_norm   = np.random.normal(loc=co_transverse, scale=transverse_spread_sigma, size=num_particles)
        p_transverse_norm = np.random.normal(scale=transverse_spread_sigma, size=num_particles)

        return [pencil, p_pencil, transverse_norm, p_transverse_norm]