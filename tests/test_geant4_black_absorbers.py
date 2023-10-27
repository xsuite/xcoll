import numpy as np
from pathlib import Path
import pytest

import xtrack as xt
import xpart as xp
import xcoll as xc
from xpart.test_helpers import flaky_assertions, retry


path = Path.cwd() / 'data_test_geant4'

npart = int(np.random.randint(1, 10) * 1e+3)

opening = 5e-3 # m

part_distribution_width = 10e-3 # m


def create_line():

    # Create a line

    # Define hor_coll_aper
    hor_coll_aper = xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-6)
    # Define hor_coll_marker
    hor_coll_marker = xt.Marker()

    # Define ver_coll_aper
    ver_coll_aper = xt.LimitRectEllipse(max_x=0.04, max_y=0.04, a_squ=0.0016, b_squ=0.0016, a_b_squ=2.56e-6)
    # Define ver_coll_marker
    ver_coll_marker = xt.Marker()

    line = xt.Line(elements=[xt.Drift(length=0.5),
                        hor_coll_aper, hor_coll_marker,
                        xt.Drift(length=1.0),
                        ver_coll_aper, ver_coll_marker,
                        xt.Drift(length=2.),
                        xt.Multipole(knl=[0, 1.], ksl=[0,0]),
                        xt.Drift(length=1.),
                        xt.Multipole(knl=[0, -1.], ksl=[0,0])],
            element_names=['drift_0', 'hor_coll_aper', 'hor_coll', 'drift_1', 'ver_coll_aper', 'ver_coll', 'drift_2', 'q1', 'drift_3', 'q2'])
    
    line.twiss_default['method'] = '4d'
    
    line.particle_ref = xp.Particles(p0c=45.6e9, q0=-1, mass0=xp.ELECTRON_MASS_EV)

    return line


def set_physical_gap(coll_manager, opening):
    
    # set horizontal collimator opening [m]
    coll_manager.line.element_dict['hor_coll'].jaw_L = opening
    coll_manager.line.element_dict['hor_coll'].jaw_R = -opening

    # set vertical collimator opening [m]
    coll_manager.line.element_dict['ver_coll'].jaw_L = opening
    coll_manager.line.element_dict['ver_coll'].jaw_R = -opening

    return coll_manager


def define_collimators():

    coll_dict = {'families':{}, 'Collimators': {}}
    coll_dict['Collimators'] = {'hor_coll' : {'angle': 0, 'length': 0.1, 'gap': 5, 'material': 'C'},
                                'ver_coll' : {'angle': 90, 'length': 0.1, 'gap': 5, 'material': 'C'}}
    
    return coll_dict


def install_black_absorbers(line, opening):

    coll_dict = define_collimators()

    # Initialize collmanager
    coll_manager = xc.CollimatorManager.from_dict(file=coll_dict, line=line, nemitt_x=0, nemitt_y=0)

    # Install collimators into line
    coll_manager.install_black_absorbers(verbose=True)

    # Build the tracker
    coll_manager.build_tracker()

    # Set the collimator openings based on the colldb,
    # or manually override with the option gaps={collname: gap}
    coll_manager.set_openings()   

    coll_manager = set_physical_gap(coll_manager, opening)

    return coll_manager 


def install_geant4_black_absorbers(line, opening):

    coll_dict = define_collimators()

    # Initialize collmanager
    coll_manager = xc.CollimatorManager.from_dict(file=coll_dict, line=line, nemitt_x=0, nemitt_y=0)

    # Install collimators into line
    coll_manager.install_geant4_collimators(verbose=True, bdsim_config_file=str(path / f'settings.gmad'))

    # Build the tracker
    coll_manager.build_tracker()

    # Set the collimator openings based on the colldb,
    # or manually override with the option gaps={collname: gap}
    coll_manager.set_openings()   

    return coll_manager 


def generate_initial_distribution(npart, part_distribution_width):

    # Generate initial distribution
    part = particles = xp.Particles(p0c=45.6e9, #eV
                                    q0=-1, mass0=xp.ELECTRON_MASS_EV,
                                    x=np.random.uniform(-part_distribution_width, part_distribution_width, npart),
                                    px=np.zeros(npart),
                                    y=np.random.uniform(-part_distribution_width, part_distribution_width, npart),
                                    py=np.zeros(npart))

    return part                                                                 


def run_black_absorbers(npart, opening, part_distribution_width):

    line = create_line()
    coll_manager = install_black_absorbers(line, opening)
    init_distrib = generate_initial_distribution(npart, part_distribution_width)
    part = init_distrib

    # Track
    coll_manager.enable_scattering()
    line.track(part, ele_start=0, ele_stop=7)
    coll_manager.disable_scattering()

    # store number of surviving particles
    n_surviving = part._num_active_particles

    return part, n_surviving, init_distrib


def run_geant4_black_absorbers(init_distrib, opening):

    line = create_line()
    coll_manager = install_geant4_black_absorbers(line, opening)
    part=init_distrib

    # Track
    coll_manager.enable_scattering()
    line.track(part, ele_start=0, ele_stop=7)
    coll_manager.disable_scattering()

    # store number of surviving particles
    n_surviving = part._num_active_particles

    return part, n_surviving


def where_losses(part, part_geant4):

    df_part = part.to_pandas()
    at_element = df_part['at_element']
    losses = at_element.value_counts()  

    df_part_geant4 = part_geant4.to_pandas()
    at_element_geant4 = df_part_geant4['at_element']
    losses_geant4 = at_element_geant4.value_counts()

    return losses, losses_geant4


def test_black_absorbers_geant4():

    part, n_surviving, init_distrib = run_black_absorbers(npart, opening, part_distribution_width)

    # run_geant4_black_absorbers takes init_distrib as input to use the same distribution used for run_black_absorbers
    part_geant4, n_surviving_geant4 = run_geant4_black_absorbers(init_distrib, opening)

    losses, losses_geant4 = where_losses(part, part_geant4)
    
    with flaky_assertions():
        assert n_surviving == n_surviving_geant4
        assert (losses.index == losses_geant4.index).all()
        assert (losses == losses_geant4).all()