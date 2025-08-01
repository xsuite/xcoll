# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from warnings import warn

import xtrack as xt
import xobjects as xo
import xpart as xp

from .beam_elements import collimator_classes, BaseCrystal


def generate_pencil_on_collimator(line, name, num_particles, *, side='+-', pencil_spread=1e-6,
                                  impact_parameter=0, sigma_z=7.61e-2, twiss=None, longitudinal=None,
                                  longitudinal_betatron_cut=None, _capacity=None, tw=None,
                                  _longitudinal_coords=None, **kwargs):
    """Generate a pencil beam on a collimator."""

    # Do some general checks
    if not line._has_valid_tracker():
        raise ValueError("Please build tracker before generating pencil distribution!")
    coll = line[name]
    if not isinstance(coll, tuple(collimator_classes)):
        raise ValueError("Need to provide a valid collimator!")
    if coll.optics is None:
        raise ValueError("Need to assign optics to collimators before generating pencil distribution!")
    num_particles = int(num_particles)
    #if _capacity is None and len(line.get_elements_of_type((Geant4Collimator, Geant4Crystal))[0]) > 0:
    #    _capacity = 2*num_particles

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

    if coll.side == 'left':
        if side == '-':
            raise ValueError("Cannot generate a pencil on the right jaw of a left-"
                             "sided collimator!")
        side = '+'
    if coll.side == 'right':
        if side == '+':
            raise ValueError("Cannot generate a pencil on the left jaw of a right-"
                             "sided collimator!")
        side = '-'

    if tw is not None:
        warn("The argument tw is deprecated. Please use twiss instead.", FutureWarning)
        if twiss is None:
            twiss = tw

    if twiss is None:
        twiss = line.twiss(reverse=False)

    # Longitudinal plane
    if _longitudinal_coords:
        zeta = _longitudinal_coords[0]
        delta = _longitudinal_coords[1]

    else:
        zeta, delta = _generate_longitudinal_dist(line, num_particles, sigma_z, longitudinal)

        # Generate 4D coordinates
        # TODO: there is some looping in the calculation here and in xpart. Can it be improved?
        if side == '+-':
            num_plus = int(num_particles/2)
            num_min  = int(num_particles - num_plus)
            zeta_plus = zeta[:num_plus] if hasattr(zeta, '__iter__') else zeta
            zeta_min  = zeta[num_plus:] if hasattr(zeta, '__iter__') else zeta
            delta_plus = delta[:num_plus] if hasattr(delta, '__iter__') else delta
            delta_min  = delta[num_plus:] if hasattr(delta, '__iter__') else delta
            if _capacity:
                _capacity_plus = int(_capacity/2)
                _capacity_min  = int(_capacity - _capacity_plus)
            else:
                _capacity_plus = None
                _capacity_min  = None
            part_plus = generate_pencil_on_collimator(line=line, name=name, num_particles=num_plus,
                                impact_parameter=impact_parameter, _capacity=_capacity_plus,
                                side='+', pencil_spread=pencil_spread, twiss=twiss,
                                _longitudinal_coords=[zeta_plus, delta_plus],
                                **kwargs)
            part_min = generate_pencil_on_collimator(line=line, name=name, num_particles=num_min,
                                impact_parameter=impact_parameter, _capacity=_capacity_min,
                                side='-', pencil_spread=pencil_spread, twiss=twiss,
                                _longitudinal_coords=[zeta_min, delta_min],
                                **kwargs)

            part = xt.Particles.merge([part_plus, part_min])
            part.start_tracking_at_element = part_plus.start_tracking_at_element
            assert part.start_tracking_at_element == part_min.start_tracking_at_element
            return part

    pencil, p_pencil, transverse_norm, p_transverse_norm, is_converging, at_element = \
                    _generate_4D_pencil_one_jaw(line, name, num_particles, plane, side,
                                                impact_parameter, pencil_spread, twiss, **kwargs)

    # Build the particles
    if plane == 'x':
        part = line.build_particles(
                x=pencil, px=p_pencil, y_norm=transverse_norm, py_norm=p_transverse_norm,
                zeta=zeta, delta=delta, nemitt_x=coll.nemitt_x, nemitt_y=coll.nemitt_y,
                at_element=at_element, _capacity=_capacity, _context=coll._buffer.context,
                **kwargs)
    else:
        part = line.build_particles(
                x_norm=transverse_norm, px_norm=p_transverse_norm, y=pencil, py=p_pencil,
                zeta=zeta, delta=delta, nemitt_x=coll.nemitt_x, nemitt_y=coll.nemitt_y,
                at_element=at_element, _capacity=_capacity, _context=coll._buffer.context,
                **kwargs)

    part._init_random_number_generator()

    if not is_converging:
        dri = xt.Drift(length=-coll.length)
        dri.track(part)
        part.start_tracking_at_element -= 1
        part.at_element -= 1

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
        twiss = line.twiss(reverse=False)

    beam_sizes = twiss.get_beam_covariance(nemitt_x=nemitt_x, nemitt_y=nemitt_y)
    beam_sizes = beam_sizes.rows[at_element:f'{at_element}>>1'][f'sigma_{plane}']
    sigma = beam_sizes[0] if match_at_front else beam_sizes[1]
    delta  = (position_mm - betatron_cut*sigma - twiss.rows[at_element][plane])
    delta /= twiss.rows[at_element][f'd{plane}']

    return delta


def _generate_4D_pencil_one_jaw(line, name, num_particles, plane, side, impact_parameter,
                                pencil_spread, twiss=None, _capacity=None, **kwargs):
    coll = line[name]
    beam_sizes = twiss.get_beam_covariance(nemitt_x=coll.nemitt_x, nemitt_y=coll.nemitt_y)

    # Is it converging or diverging?
    # TODO: dispersion might change this
    # TODO: skew collimators
    tolerance_tilt = 1e-12 # 0.1 urad tolerance on jaw tilt  =>  we prioritise converging
    divergence = coll.divergence
    if side == '+':
        if isinstance(coll, BaseCrystal):
            # A pencil on the crystal should always be upstream
            is_converging = True
            pencil_pos = coll.jaw_U + impact_parameter
        else:
            betatron_angle = coll.gap_L * divergence
            is_converging = coll.tilt_L + tolerance_tilt >= betatron_angle
            print(f"Left jaw of collimator {name} is {'con' if is_converging else 'di'}verging.")
            if is_converging:
                pencil_pos = coll.jaw_LU + impact_parameter
            else:
                pencil_pos = coll.jaw_LD + impact_parameter
    elif side == '-':
        if isinstance(coll, BaseCrystal):
            # A pencil on the crystal should always be upstream
            is_converging = True
            pencil_pos = coll.jaw_U - impact_parameter
        else:
            betatron_angle = coll.gap_R * divergence
            is_converging = coll.tilt_R - tolerance_tilt <= betatron_angle
            print(f"Right jaw of collimator {name} is {'con' if is_converging else 'di'}verging.")
            if is_converging:
                pencil_pos = coll.jaw_RU - impact_parameter
            else:
                pencil_pos = coll.jaw_RD - impact_parameter
    else:
        raise ValueError(f"Sinde {side} not supported in  _generate_4D_pencil_one_jaw!")

    if is_converging:
        # pencil at front of jaw
        sigma = beam_sizes.rows[name:f'{name}>>1'][f'sigma_{plane}'][0]
        tw_at_s = twiss.rows[name]
        at_element = name
    else:
        # pencil at back of jaw
        sigma = beam_sizes.rows[name:f'{name}>>1'][f'sigma_{plane}'][1]
        tw_at_s = twiss.rows[f'{name}>>1']
        at_element = line.element_names[line.element_names.index(name)+1]

    dr_sigmas = pencil_spread/sigma

    # Collimator plane: generate pencil distribution
    pencil, p_pencil = xp.generate_2D_pencil_with_absolute_cut(
                        num_particles, plane=plane, absolute_cut=pencil_pos, line=line,
                        dr_sigmas=dr_sigmas, nemitt_x=coll.nemitt_x, nemitt_y=coll.nemitt_y,
                        at_element=at_element, side=side, twiss=tw_at_s, **kwargs
                       )

    # Other plane: generate gaussian distribution in normalized coordinates
    transverse_norm   = np.random.normal(size=num_particles)
    p_transverse_norm = np.random.normal(size=num_particles)

    return pencil, p_pencil, transverse_norm, p_transverse_norm, is_converging, at_element


def _generate_longitudinal_dist(line, num_particles, sigma_z, longitudinal):
    # TODO: make this more general, make this better
    if longitudinal is None:
        return 0, 0
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
        return xp.generate_longitudinal_coordinates(
                num_particles=num_particles, distribution='gaussian', sigma_z=sigma_z, line=line
        )
    elif not hasattr(longitudinal, '__iter__'):
        raise ValueError
    elif len(longitudinal) != 2:
        raise ValueError
    elif isinstance(longitudinal, str):
        raise ValueError
    elif isinstance(longitudinal, dict):
        return longitudinal['zeta'], longitudinal['delta']
    else:
        return longitudinal[0], longitudinal[1]
