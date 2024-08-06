# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np

import xtrack as xt
import xobjects as xo
import xpart as xp

from .beam_elements import collimator_classes


def generate_pencil_on_collimator(line, name, num_particles, *, side='+-', pencil_spread=1e-6,
                                  impact_parameter=0, sigma_z=7.61e-2, tw=None, longitudinal=None,
                                  longitudinal_betatron_cut=None, _capacity=None):
    """
    Generate a pencil beam on a collimator.
    """

    if not line._has_valid_tracker():
        raise Exception("Please build tracker before generating pencil distribution!")

    coll = line[name]

    if not isinstance(coll, tuple(collimator_classes)):
        raise ValueError("Need to provide a valid collimator!")

    if coll.optics is None:
        raise Exception("Need to assign optics to collimators before generating pencil distribution!")

    num_particles = int(num_particles)

    if coll.side == 'left':
        side = '+'
    if coll.side == 'right':
        side = '-'

    # Define the plane
    angle = coll.angle
    if abs(np.mod(angle-90,180)-90) < 1e-6:
        plane = 'x'
        transv_plane = 'y'
    elif abs(np.mod(angle,180)-90) < 1e-6:
        plane = 'y'
        transv_plane = 'x'
    else:
        raise NotImplementedError("Pencil beam on a skew collimator not yet supported!")

    if tw is None:
        tw = line.twiss()    # TODO: can we do this smarter by caching?

    # Is it converging or diverging?    # TODO: This might change with a tilt!!!!!!
    s_front = line.get_s_position(name)
    s_back  = s_front + coll.length
    is_converging  = tw[f'alf{plane}', name] > 0
    print(f"Collimator {name} is {'con' if is_converging else 'di'}verging.")

    beam_sizes = tw.get_beam_covariance(nemitt_x=coll.nemitt_x, nemitt_y=coll.nemitt_y)
    if is_converging:
        # pencil at front of jaw
        match_at_s = s_front
        sigma = beam_sizes.rows[name:f'{name}%%1'][f'sigma_{plane}'][0]
        sigma_transv = beam_sizes.rows[name:f'{name}%%1'][f'sigma_{transv_plane}'][0]
    else:
        # pencil at back of jaw
        match_at_s = s_back
        sigma = beam_sizes.rows[name:f'{name}%%1'][f'sigma_{plane}'][1]
        sigma_transv = beam_sizes.rows[name:f'{name}%%1'][f'sigma_{transv_plane}'][1]
    dr_sigmas = pencil_spread/sigma

    # Generate 4D coordinates
    # TODO: there is some looping in the calculation here and in xpart. Can it be improved?
    if side == '+-':
        num_plus = int(num_particles/2)
        num_min  = int(num_particles - num_plus)
        coords_plus = _generate_4D_pencil_one_jaw(line, name, num_plus, plane, '+', impact_parameter, dr_sigmas, match_at_s, is_converging)
        coords_min  = _generate_4D_pencil_one_jaw(line, name, num_min,  plane, '-', impact_parameter, dr_sigmas, match_at_s, is_converging)
        coords      = [ [*c_plus, *c_min] for c_plus, c_min in zip(coords_plus, coords_min)]
    else:
        coords      = _generate_4D_pencil_one_jaw(line, name, num_particles, plane, side, impact_parameter, dr_sigmas, match_at_s, is_converging)
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
        raise NotImplementedError
        # if longitudinal_betatron_cut is None:
        #     cut = 0
        # else:
        #     cut = np.random.uniform(-longitudinal_betatron_cut, longitudinal_betatron_cut,
        #                             num_particles)
        # delta = generate_delta_from_dispersion(line, name, plane=plane, position_mm=pencil,
        #                                        nemitt_x=nemitt_x, nemitt_y=nemitt_y, twiss=tw,
        #                                        betatron_cut=cut, match_at_front=is_converging)
        # zeta  = 0
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
                zeta=zeta, delta=delta, nemitt_x=coll.nemitt_x, nemitt_y=coll.nemitt_y,
                line=line, at_element=name, match_at_s=match_at_s,
                _context=coll._buffer.context, _capacity=_capacity
        )
    else:
        part = xp.build_particles(
                x_norm=transverse_norm, px_norm=p_transverse_norm, y=pencil, py=p_pencil, 
                zeta=zeta, delta=delta, nemitt_x=coll.nemitt_x, nemitt_y=coll.nemitt_y,
                line=line, at_element=name, match_at_s=match_at_s,
                _context=coll._buffer.context, _capacity=_capacity
        )

    part._init_random_number_generator()

    return part


def generate_delta_from_dispersion(line, at_element, *, plane, position_mm, nemitt_x, nemitt_y,
                                   twiss=None, betatron_cut=0, match_at_front=True):
    if line.tracker is None:
        raise ValueError("Need to build tracker first!")
    if not hasattr(betatron_cut, '__iter__'):
        if hasattr(position_mm, '__iter__'):
            betatron_cut = np.full_like(position_mm, betatron_cut)
    elif not hasattr(position_mm, '__iter__'):
        position_mm = np.full_like(betatron_cut, position_mm)
    elif len(position_mm) != len(betatron_cut):
        raise ValueError
    if plane not in ['x', 'y']:
        raise ValueError("The variable 'plane' needs to be either 'x' or 'y'!")

    if twiss is None:
        twiss = line.twiss()

    beam_sizes = twiss.get_beam_covariance(nemitt_x=nemitt_x, nemitt_y=nemitt_y)
    beam_sizes = beam_sizes.rows[at_element:f'{at_element}%%1'][f'sigma_{plane}']
    sigma = beam_sizes[0] if match_at_front else beam_sizes[1]
    delta  = (position_mm - betatron_cut*sigma - twiss.rows[at_element][plane])
    delta /= twiss.rows[at_element][f'd{plane}']

    return delta


def _generate_4D_pencil_one_jaw(line, name, num_particles, plane, side, impact_parameter,
                                dr_sigmas, match_at_s, is_converging):
    coll = line[name]

    if side == '+':
        if is_converging:
            pencil_pos = coll.jaw_LU + impact_parameter
        else:
            pencil_pos = coll.jaw_LD + impact_parameter
    elif side == '-':
        if is_converging:
            pencil_pos = coll.jaw_RU - impact_parameter
        else:
            pencil_pos = coll.jaw_RD - impact_parameter

    # Collimator plane: generate pencil distribution
    pencil, p_pencil = xp.generate_2D_pencil_with_absolute_cut(
                        num_particles, plane=plane, absolute_cut=pencil_pos, line=line,
                        dr_sigmas=dr_sigmas, nemitt_x=coll.nemitt_x, nemitt_y=coll.nemitt_y,
                        at_element=name, side=side,match_at_s=match_at_s
                       )

    # Other plane: generate gaussian distribution in normalized coordinates
    transverse_norm   = np.random.normal(size=num_particles)
    p_transverse_norm = np.random.normal(size=num_particles)

    return [pencil, p_pencil, transverse_norm, p_transverse_norm]

