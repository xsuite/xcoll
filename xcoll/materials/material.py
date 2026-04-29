# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
import scipy.constants as sc
from scipy.interpolate import CubicSpline
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional
from pathlib import Path
import xobjects as xo

from ..general import _pkg_root
from .parameters import (_approximate_radiation_length, _default_excitation_energies,
                         _combine_radiation_lengths, _combine_excitation_energies,
                         _average_Z_over_A, _effective_Z2)
from ..compare import deep_equal
from .glauber_gribov import SPECIES, load_all_splines, make_nucleon_cs
from .isotopes import ISOTOPES
from .build_isotopes import SYMBOL_TO_Z
_materials_context = xo.ContextCpu()

# CLASS INHERITANCE DOES NOT WORK WITH DYNAMIC XOFIELDS
# Because then the memory offsets between parent and child might no longer be the same,
# as xobjects enforces a strict order: static fields first, and then dynamic fields.
# See struct.py, in __new__ of MetaStruct

Z_TO_SYMBOL = {v: k for k, v in SYMBOL_TO_Z.items()}

N_CS_POINTS = 1000
@dataclass(frozen=True)
class IsotopeData:
    """
    One isotope of an element.
 
    Attributes
    ----------
    mass_number  : integer A
    atomic_mass  : precise atomic mass in unified atomic mass units [u]
    abundance    : natural isotopic abundance (fraction, sums to 1 per element).
                   None for isotopes with no stable natural occurrence.
    """
    mass_number: int
    atomic_mass: float
    abundance:   Optional[float]

class Material(xo.HybridClass):
    _xofields = {
        # Essential fields
        '_density':                 xo.Float64,     # [g/cm^3], can be auto-calculated for compounds
        # Auto-calculated fields
        '_ZA_mean':                 xo.Float64,     # [mol of electrons per gram]
        '_Z2_eff':                  xo.Float64,     # Effective Z for Rutherford scattering
        '_atoms_per_volume':        xo.Float64,     # [atoms/m^3]
        '_num_nucleons_eff':        xo.Float64,     # Effective number of nucleons for nuclear interactions

        # Cross sections
        '_cs_knots_pp':          xo.Float64[N_CS_POINTS],       # log(sqrt_s) knot positions
        '_n_points_pp':             xo.Float64,                      # Number of points in the GG arrays
        '_cs_sqrt_s_pp':         xo.Float64[N_CS_POINTS],       # Minimum sqrt(s) for the spline
        '_cs_log_sqrt_s_min_pp': xo.Float64,                    # Minimum log(sqrt(s)) for the spline
        '_cs_log_step_pp':       xo.Float64,                    # Step size for the spline
        '_cs_knots_kmin':          xo.Float64[N_CS_POINTS],       # log(sqrt_s) knot positions
        '_n_points_kmin':          xo.Float64,                      # Number of points in the GG arrays
        '_cs_sqrt_s_kmin':         xo.Float64[N_CS_POINTS],       # Minimum sqrt(s) for the spline
        '_cs_log_sqrt_s_min_kmin': xo.Float64,                    # Minimum log(sqrt(s)) for the spline
        '_cs_log_step_kmin':       xo.Float64,                    # Step size for the spline
        '_cs_knots_kplus':          xo.Float64[N_CS_POINTS],       # log(sqrt_s) knot positions
        '_n_points_kplus':          xo.Float64,                      # Number of points in the GG arrays
        '_cs_sqrt_s_kplus':         xo.Float64[N_CS_POINTS],       # Minimum sqrt(s) for the spline
        '_cs_log_sqrt_s_min_kplus': xo.Float64,                    # Minimum log(sqrt(s)) for the spline
        '_cs_log_step_kplus':       xo.Float64,                    # Step size for the spline
        '_cs_knots_pimin':          xo.Float64[N_CS_POINTS],       # log(sqrt_s) knot positions
        '_n_points_pimin':          xo.Float64,                      # Number of points in the GG arrays
        '_cs_sqrt_s_pimin':         xo.Float64[N_CS_POINTS],       # Minimum sqrt(s) for the spline
        '_cs_log_sqrt_s_min_pimin': xo.Float64,                    # Minimum log(sqrt(s)) for the spline
        '_cs_log_step_pimin':       xo.Float64,                    # Step size for the spline
        '_cs_knots_piplus':          xo.Float64[N_CS_POINTS],       # log(sqrt_s) knot positions
        '_n_points_piplus':          xo.Float64,                      # Number of points in the GG arrays
        '_cs_sqrt_s_piplus':         xo.Float64[N_CS_POINTS],       # Minimum sqrt(s) for the spline
        '_cs_log_sqrt_s_min_piplus': xo.Float64,                    # Minimum log(sqrt(s)) for the spline
        '_cs_log_step_piplus':       xo.Float64,                    # Step size for the spline

        '_nuclear_slope':      xo.Float64,                    # Nuclear slope parameter
        '_nuclear_slope_pion': xo.Float64,                    # Nuclear slope parameter for pions
        # # GG raw values (for diagnostics)
        # '_cs_tot_pp_GG':        xo.Float64[N_CS_POINTS],
        # '_cs_el_pp_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_prod_pp_GG':       xo.Float64[N_CS_POINTS],
        # '_cs_sd_pp_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_el_nucleon_pp_GG':    xo.Float64[N_CS_POINTS],

        # '_cs_tot_kmin_GG':        xo.Float64[N_CS_POINTS],
        # '_cs_el_kmin_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_prod_kmin_GG':       xo.Float64[N_CS_POINTS],
        # '_cs_sd_kmin_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_el_nucleon_kmin_GG':    xo.Float64[N_CS_POINTS],

        # '_cs_tot_kplus_GG':        xo.Float64[N_CS_POINTS],
        # '_cs_el_kplus_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_prod_kplus_GG':       xo.Float64[N_CS_POINTS],
        # '_cs_sd_kplus_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_el_nucleon_kplus_GG':    xo.Float64[N_CS_POINTS],

        # '_cs_tot_pimin_GG':        xo.Float64[N_CS_POINTS],
        # '_cs_el_pimin_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_prod_pimin_GG':       xo.Float64[N_CS_POINTS],
        # '_cs_sd_pimin_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_el_nucleon_pimin_GG':    xo.Float64[N_CS_POINTS],

        # '_cs_tot_piplus_GG':        xo.Float64[N_CS_POINTS],
        # '_cs_el_piplus_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_prod_piplus_GG':       xo.Float64[N_CS_POINTS],
        # '_cs_sd_piplus_GG':         xo.Float64[N_CS_POINTS],
        # '_cs_el_nucleon_piplus_GG':    xo.Float64[N_CS_POINTS],
        # Spline coefficients — shape (N_CS_POINTS - 1,) per coefficient
        '_cs_tot_pp_a':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_pp_b':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_pp_c':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_pp_d':      xo.Float64[N_CS_POINTS-1],
        '_cs_el_pp_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_pp_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_pp_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_pp_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_prod_pp_a':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_pp_b':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_pp_c':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_pp_d':     xo.Float64[N_CS_POINTS-1],
        '_cs_sd_pp_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_pp_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_pp_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_pp_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_pp_a':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_pp_b':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_pp_c':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_pp_d':  xo.Float64[N_CS_POINTS-1],

        '_cs_tot_kmin_a':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_kmin_b':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_kmin_c':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_kmin_d':      xo.Float64[N_CS_POINTS-1],
        '_cs_el_kmin_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_kmin_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_kmin_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_kmin_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_prod_kmin_a':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_kmin_b':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_kmin_c':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_kmin_d':     xo.Float64[N_CS_POINTS-1],
        '_cs_sd_kmin_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_kmin_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_kmin_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_kmin_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_kmin_a':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_kmin_b':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_kmin_c':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_kmin_d':  xo.Float64[N_CS_POINTS-1],

        '_cs_tot_kplus_a':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_kplus_b':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_kplus_c':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_kplus_d':      xo.Float64[N_CS_POINTS-1],
        '_cs_el_kplus_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_kplus_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_kplus_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_kplus_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_prod_kplus_a':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_kplus_b':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_kplus_c':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_kplus_d':     xo.Float64[N_CS_POINTS-1],
        '_cs_sd_kplus_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_kplus_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_kplus_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_kplus_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_kplus_a':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_kplus_b':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_kplus_c':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_kplus_d':  xo.Float64[N_CS_POINTS-1],

        '_cs_tot_pimin_a':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_pimin_b':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_pimin_c':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_pimin_d':      xo.Float64[N_CS_POINTS-1],
        '_cs_el_pimin_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_pimin_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_pimin_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_pimin_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_prod_pimin_a':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_pimin_b':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_pimin_c':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_pimin_d':     xo.Float64[N_CS_POINTS-1],
        '_cs_sd_pimin_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_pimin_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_pimin_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_pimin_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_pimin_a':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_pimin_b':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_pimin_c':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_pimin_d':  xo.Float64[N_CS_POINTS-1],

        '_cs_tot_piplus_a':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_piplus_b':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_piplus_c':      xo.Float64[N_CS_POINTS-1],
        '_cs_tot_piplus_d':      xo.Float64[N_CS_POINTS-1],
        '_cs_el_piplus_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_piplus_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_piplus_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_piplus_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_prod_piplus_a':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_piplus_b':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_piplus_c':     xo.Float64[N_CS_POINTS-1],
        '_cs_prod_piplus_d':     xo.Float64[N_CS_POINTS-1],
        '_cs_sd_piplus_a':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_piplus_b':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_piplus_c':       xo.Float64[N_CS_POINTS-1],
        '_cs_sd_piplus_d':       xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_piplus_a':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_piplus_b':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_piplus_c':  xo.Float64[N_CS_POINTS-1],
        '_cs_el_nucleon_piplus_d':  xo.Float64[N_CS_POINTS-1],
        # Auto-calculated fields but can be provided for more precision
        '_radiation_length':        xo.Float64,     # [m]
        '_excitation_energy':       xo.Float64,     # [eV]
        # Optional fields (needed for full Everest support)
        '_nuclear_radius':          xo.Float64,     # emr
        '_nuclear_elastic_slope':   xo.Float64,     # [g/cm^2]    ~ 14.1 A^0.65
                # slope from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.21.3010
        '_cross_section':           xo.Float64[6],  # [barn]     ~ these combine linear in atomic fractions
                # Index 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
                #       4:Single Diffractive pp or pn, 5:Coulomb for t above mcs
        '_hcut':                    xo.Float64,     # Cut in Rutherford distribution
        # Optional crystal material fields (needed for full Everest crystal support)
        '_crystal_plane_distance':   xo.Float64,    # ai  [mm]
        '_crystal_potential':        xo.Float64,    # eum  [eV]
        '_nuclear_collision_length': xo.Float64,    # collnt [m]
        '_eta':                      xo.Float64
    }

    _depends_on = [xo.Float64] # A HybridClass needs something to depend on, otherwise the class is added twice in the cdefs during compilation

    _skip_in_to_dict  = ['_ZA_mean', '_Z2_eff', '_density', '_radiation_length',
                         '_excitation_energy', '_atoms_per_volume',
                         '_num_nucleons_eff', '_nuclear_radius',
                         '_nuclear_elastic_slope', '_cross_section',
                         '_crystal_plane_distance', '_crystal_potential',
                         '_nuclear_collision_length', '_eta',
                         '_hcut', '_components']
    _store_in_to_dict = ['Z', 'A', 'components', 'n_atoms', 'mass_fractions',
                         'density', 'radiation_length', 'excitation_energy',
                         'nuclear_radius', 'nuclear_elastic_slope',
                         'crystal_plane_distance', 'crystal_potential',
                         'nuclear_collision_length', 'eta',
                         'cross_section', 'hcut', 'state', 'temperature',
                         'pressure', 'info', 'name', 'short_name',
                         'geant4_name', 'fluka_name','isotopes']

    _kernels = {'evaluate_glauber_spline': xo.Kernel(
                                c_name='MaterialData_evaluate_glauber_spline',
                                args=[xo.Arg(xo.ThisClass, name="material"),
                                      xo.Arg(xo.Float64, name="sqrt_s"),
                                      xo.Arg(xo.Int8, name="key"),
                                      xo.Arg(xo.Int8, name="particle_id")],
                                ret=xo.Arg(xo.Float64, name="evaluate_glauber_spline")),
                }
    _needs_compilation = True
    _extra_c_sources = [
        _pkg_root / 'materials' / 'glauber_splines.h'
    ]
    # ======================
    # === Initialisation ===
    # ======================
    def __sh__(self):
        return id(self)

    def __eq__(self, other):
        return self is other    
    def __init__(self, **kwargs):
        if '_xobject' in kwargs and kwargs['_xobject'] is not None:
            super().__init__(**kwargs)
            return

        # Create xobject with all invalid values (-1)
        xokwargs = kwargs.pop('_xokwargs', {})
        for kk in ('_ZA_mean', '_Z2_eff', '_nuclear_slope', '_nuclear_slope_pion', '_radiation_length', '_excitation_energy',
                   '_atoms_per_volume', '_num_nucleons_eff', '_density', '_nuclear_radius',
                   '_nuclear_elastic_slope', '_hcut', '_crystal_plane_distance',
                   '_crystal_potential', '_eta', '_nuclear_collision_length',
                   '_cs_log_sqrt_s_min_pp', '_cs_log_step_pp', '_n_points_pp',
                   '_cs_log_sqrt_s_min_kmin', '_cs_log_step_kmin', '_n_points_kmin',
                   '_cs_log_sqrt_s_min_kplus', '_cs_log_step_kplus', '_n_points_kplus',
                   '_cs_log_sqrt_s_min_pimin', '_cs_log_step_pimin', '_n_points_pimin',
                   '_cs_log_sqrt_s_min_piplus', '_cs_log_step_piplus', '_n_points_piplus'):
            xokwargs[kk] = kwargs.pop(kk, -1.)

        for kk in (
                #    '_cs_tot_pp_GG', '_cs_el_pp_GG', '_cs_prod_pp_GG', '_cs_sd_pp_GG', '_cs_el_nucleon_pp_GG', 
                #    '_cs_tot_kmin_GG', '_cs_el_kmin_GG', '_cs_prod_kmin_GG', '_cs_sd_kmin_GG', '_cs_el_nucleon_kmin_GG',
                #    '_cs_tot_kplus_GG', '_cs_el_kplus_GG', '_cs_prod_kplus_GG', '_cs_sd_kplus_GG', '_cs_el_nucleon_kplus_GG',
                #    '_cs_tot_pimin_GG', '_cs_el_pimin_GG', '_cs_prod_pimin_GG', '_cs_sd_pimin_GG', '_cs_el_nucleon_pimin_GG',
                #    '_cs_tot_piplus_GG', '_cs_el_piplus_GG', '_cs_prod_piplus_GG', '_cs_sd_piplus_GG', '_cs_el_nucleon_piplus_GG',
        
                   '_cs_tot_pp', '_cs_el_pp','_cs_inel_pp', '_cs_tot_pn', 'cs_tot_kmin_p', 'cs_el_kmin_p', 'cs_inel_kmin_p', 
                   'cs_tot_kplus_p', 'cs_el_kplus_p', 'cs_inel_kplus_p','cs_tot_pimin_p', 'cs_el_pimin_p', 'cs_inel_pimin_p',
                   'cs_tot_piplus_p', 'cs_el_piplus_p', 'cs_inel_piplus_p', 
    
                   '_cs_sqrt_s_pp', '_cs_sqrt_s_kmin', '_cs_sqrt_s_kplus', '_cs_sqrt_s_pimin', '_cs_sqrt_s_piplus',
                   '_cs_knots_pp', '_cs_knots_kmin', '_cs_knots_kplus', '_cs_knots_pimin', '_cs_knots_piplus'):
            xokwargs[kk] = kwargs.pop(kk, [-1.] * N_CS_POINTS)

        for kk in ('_cs_tot_pp_a',  '_cs_tot_pp_b',  '_cs_tot_pp_c',  '_cs_tot_pp_d',
                   '_cs_el_pp_a',   '_cs_el_pp_b',   '_cs_el_pp_c',   '_cs_el_pp_d',
                   '_cs_prod_pp_a', '_cs_prod_pp_b', '_cs_prod_pp_c', '_cs_prod_pp_d',
                   '_cs_sd_pp_a',   '_cs_sd_pp_b',   '_cs_sd_pp_c',   '_cs_sd_pp_d',
                   '_cs_el_nucleon_pp_a',  '_cs_el_nucleon_pp_b',  '_cs_el_nucleon_pp_c',  '_cs_el_nucleon_pp_d',
                   '_cs_tot_kmin_a',  '_cs_tot_kmin_b',  '_cs_tot_kmin_c',  '_cs_tot_kmin_d',
                   '_cs_el_kmin_a',   '_cs_el_kmin_b',   '_cs_el_kmin_c',   '_cs_el_kmin_d',
                   '_cs_prod_kmin_a', '_cs_prod_kmin_b', '_cs_prod_kmin_c', '_cs_prod_kmin_d',
                   '_cs_sd_kmin_a',   '_cs_sd_kmin_b',   '_cs_sd_kmin_c',   '_cs_sd_kmin_d',
                   '_cs_el_nucleon_kmin_a',  '_cs_el_nucleon_kmin_b',  '_cs_el_nucleon_kmin_c',  '_cs_el_nucleon_kmin_d',
                   '_cs_tot_kplus_a',  '_cs_tot_kplus_b',  '_cs_tot_kplus_c',  '_cs_tot_kplus_d',
                   '_cs_el_kplus_a',   '_cs_el_kplus_b',   '_cs_el_kplus_c',   '_cs_el_kplus_d',
                   '_cs_prod_kplus_a', '_cs_prod_kplus_b', '_cs_prod_kplus_c', '_cs_prod_kplus_d',
                   '_cs_sd_kplus_a',   '_cs_sd_kplus_b',   '_cs_sd_kplus_c',   '_cs_sd_kplus_d',
                   '_cs_el_nucleon_kplus_a',  '_cs_el_nucleon_kplus_b',  '_cs_el_nucleon_kplus_c',  '_cs_el_nucleon_kplus_d',
                   '_cs_tot_pimin_a',  '_cs_tot_pimin_b',  '_cs_tot_pimin_c',  '_cs_tot_pimin_d',
                   '_cs_el_pimin_a',   '_cs_el_pimin_b',   '_cs_el_pimin_c',   '_cs_el_pimin_d',
                   '_cs_prod_pimin_a', '_cs_prod_pimin_b', '_cs_prod_pimin_c', '_cs_prod_pimin_d',
                   '_cs_sd_pimin_a',   '_cs_sd_pimin_b',   '_cs_sd_pimin_c',   '_cs_sd_pimin_d',
                   '_cs_el_nucleon_pimin_a',  '_cs_el_nucleon_pimin_b',  '_cs_el_nucleon_pimin_c',  '_cs_el_nucleon_pimin_d',
                   '_cs_tot_piplus_a',  '_cs_tot_piplus_b',  '_cs_tot_piplus_c',  '_cs_tot_piplus_d',
                   '_cs_el_piplus_a',   '_cs_el_piplus_b',   '_cs_el_piplus_c',   '_cs_el_piplus_d',
                   '_cs_prod_piplus_a', '_cs_prod_piplus_b', '_cs_prod_piplus_c', '_cs_prod_piplus_d',
                   '_cs_sd_piplus_a',   '_cs_sd_piplus_b',   '_cs_sd_piplus_c',   '_cs_sd_piplus_d',
                   '_cs_el_nucleon_piplus_a',  '_cs_el_nucleon_piplus_b',  '_cs_el_nucleon_piplus_c',  '_cs_el_nucleon_piplus_d'):
            xokwargs[kk] = kwargs.pop(kk, [-1.] * (N_CS_POINTS - 1))

        xokwargs['_cross_section'] = kwargs.pop('_cross_section', [-1., -1., -1., -1., -1., -1.])
        xokwargs['_context'] = kwargs.pop('_context', _materials_context)  # This is needed to get all materials in the same buffer (otherwise Xtrack tests fail)
        xokwargs['__class__'] = kwargs.pop('__class__', Material)
        xokwargs['_buffer'] = kwargs.pop('_buffer', None)
        xokwargs['_offset'] = kwargs.pop('_offset', None)
        xokwargs['_kwargs_name_check'] = kwargs.pop('_kwargs_name_check', None)
        super().__init__(**xokwargs)

        # Set python-side defaults
        self._Z = None
        self._A = None
        self._components = None
        self._n_atoms = None
        self._mass_fractions = None
        self._radiation_length_set_manually = False
        self._excitation_energy_set_manually = False
        self._state = None
        self._temperature = None
        self._pressure = None
        self._info = None
        self._name = None
        self._short_name = None
        self._fluka_name = None
        self._geant4_name = None
        self._out_of_sync = False
        self._frozen = False  # Pre-defined materials will be frozen at package import
        self._generated_geant4_code = None
        self._generated_fluka_code = None
        self._isotopes = None
        self.isotopes = kwargs.pop('isotopes', None)

        # For the mandatory fields, decide how to initialise (elemental or compound)
        if ('Z' in kwargs or 'A' in kwargs) and ('components' in kwargs \
        or 'n_atoms' in kwargs or 'mass_fractions' in kwargs
        or 'volume_fractions' in kwargs or 'molar_fractions' in kwargs \
        or 'atomic_fractions' in kwargs):
            raise ValueError("Invalid material definition! Use either `Z` "
                    "and `A` for elemental materials, or `components` and "
                    "`n_atoms`, `mass_fractions`, `volume_fractions`, "
                    "`molar_fractions`, or `atomic_fractions` for compound "
                    "materials.")

        if 'Z' in kwargs or 'A' in kwargs:
            self._init_element(kwargs)

        elif 'components' in kwargs or 'n_atoms' in kwargs \
        or 'mass_fractions' in kwargs or 'volume_fractions' in kwargs \
        or 'molar_fractions' in kwargs or 'atomic_fractions' in kwargs:
            self._init_compound(kwargs)

        else:
            raise ValueError("Invalid material definition! Use either `Z` "
                    "and `A` for elemental materials, or `components` and "
                    "`n_atoms`, `mass_fractions`, `volume_fractions`, "
                    "`molar_fractions`, or `atomic_fractions` for compound "
                    "materials.")

        # Required properties
        if 'density' not in kwargs:
            raise ValueError('density must be provided for Material')
        self.density = kwargs.pop('density')
        self.radiation_length = kwargs.pop('radiation_length', None)   # Can be provided for more precision
        self.excitation_energy = kwargs.pop('excitation_energy', None) # Can be provided for more precision
        self.update_vars()

        # Assign optional properties
        for kk in ['nuclear_radius', 'nuclear_elastic_slope', 'cross_section', 'hcut',
                   'crystal_plane_distance', 'crystal_potential', 'eta',
                   'nuclear_collision_length']:
            if kk in kwargs:
                setattr(self, kk, kwargs.pop(kk))
        for kk in ['state', 'temperature', 'pressure', 'info']:
            if kk in kwargs:
                setattr(self, kk, kwargs.pop(kk))
        for kk in ['name', 'short_name', 'fluka_name', 'geant4_name']:
            if kk in kwargs:
                setattr(self, kk, kwargs.pop(kk))

        # Check for unused kwargs
        if len(kwargs) > 0:
            raise ValueError(f"Unknown keyword arguments in Material "
                             f"initialisation: {list(kwargs.keys())}")


    def _init_element(self, kwargs):
        self.Z = kwargs.pop('Z', None)
        self.A = kwargs.pop('A', None)
        if self.Z is None:
            raise ValueError('Z must be provided for an elemental Material')
        if self.A is None:
            raise ValueError('A must be provided for an elemental Material')
        self._get_material_cross_sections()

    def _init_compound(self, kwargs):
        self._resolve_components(kwargs)
        self._resolve_elements(kwargs)
        if any([el.name is None for el in self.components]):
            raise ValueError('All components must have a name')
        for el in self.components:
            if el._out_of_sync:
                raise ValueError(f"Component material {el.name} is out of sync "
                                 f"with database. Cannot use it to build a "
                                 f"compound.")
        self._get_material_cross_sections()

    def _resolve_components(self, kwargs):
        components = kwargs.pop('components', None)
        if components is None:
            components = kwargs.pop('_components', None)
            if components is None:
                raise ValueError('Variable `components` must be provided')
        if len(components) < 2:
            raise ValueError('Material must have at least two components')
        from xcoll.materials.database import db as mdb
        components = [mdb[el] if isinstance(el, str) else el for el in components]
        components = [Material.from_dict(el) if isinstance(el, dict) else el for el in components]
        if any([not isinstance(el, Material) for el in components]):
            raise ValueError('All components must be of type Material')
        self._components = components

    def _resolve_elements(self, kwargs):
        n_atoms = kwargs.pop('n_atoms', None)
        mass_fractions = kwargs.pop('mass_fractions', None)
        volume_fractions = kwargs.pop('volume_fractions', None)
        molar_fractions = kwargs.pop('molar_fractions', None)
        atomic_fractions = kwargs.pop('atomic_fractions', None)
        if n_atoms is not None:
            # Composition is defined by number of atoms of each element
            if mass_fractions is not None:
                raise ValueError("Cannot provide both `mass_fractions` "
                                 "and `n_atoms`")
            if volume_fractions is not None:
                raise ValueError("Cannot provide both `volume_fractions` "
                                 "and `n_atoms`")
            if molar_fractions is not None:
                raise ValueError("Cannot provide both `molar_fractions` "
                                 "and `n_atoms`")
            if atomic_fractions is not None:
                raise ValueError("Cannot provide both `atomic_fractions` "
                                 "and `n_atoms`")
            self._resolve_n_atoms(self.components, n_atoms)
        elif volume_fractions is not None:
            # Composition is defined by volume fractions of each element
            if mass_fractions is not None:
                raise ValueError("Cannot provide both `mass_fractions` "
                                 "and `volume_fractions`")
            if molar_fractions is not None:
                raise ValueError("Cannot provide both `molar_fractions` "
                                 "and `volume_fractions`")
            if atomic_fractions is not None:
                raise ValueError("Cannot provide both `atomic_fractions` "
                                 "and `volume_fractions`")
            self._resolve_volume_fractions(self.components, volume_fractions)
        elif molar_fractions is not None:
            # Composition is defined by molar fractions of each element
            if mass_fractions is not None:
                raise ValueError("Cannot provide both `mass_fractions` "
                                 "and `molar_fractions`")
            if atomic_fractions is not None:
                raise ValueError("Cannot provide both `atomic_fractions` "
                                 "and `molar_fractions`")
            self._resolve_molar_fractions(self.components, molar_fractions)
        elif atomic_fractions is not None:
            # Composition is defined by atomic fractions of each element
            if mass_fractions is not None:
                raise ValueError("Cannot provide both `mass_fractions` "
                                 "and `atomic_fractions`")
            self._resolve_molar_fractions(self.components, atomic_fractions)
        elif mass_fractions is not None:
            # Composition is defined by mass fractions of each element
            self._resolve_mass_fractions(self.components, mass_fractions)
        else:
            raise ValueError("One of `n_atoms`, `mass_fractions`, "
                             "`volume_fractions`, or `molar_fractions` "
                             "must be provided")

    def _resolve_n_atoms(self, components, n_atoms):
        if any([nn <= 0 for nn in n_atoms]):
            raise ValueError('All n_atoms must be strictly positive')
        while any([el.components is not None for el in components]):
            this_components = []
            this_n_atoms = []
            for el, nn in zip(components, n_atoms):
                if el.components is not None:
                    if not el.n_atoms:
                        raise ValueError("When defining a Material by `n_atoms`, "
                                            "any nested Material must have `n_atoms` "
                                            "defined as well.")
                    this_components.extend(el.components)
                    this_n_atoms.extend([n * nn for n in el.n_atoms])
                else:
                    this_components.append(el)
                    this_n_atoms.append(nn)
            components = this_components
            n_atoms = this_n_atoms
        if len(components) != len(n_atoms):
            raise ValueError('Variables `components` and `n_atoms` must have the same length')
        # Sum duplicates
        agg = defaultdict(float)
        for el, nn in zip(components, n_atoms):
            agg[el] += nn
        components, n_atoms = zip(*agg.items())
        self._components = np.array(components)
        self._n_atoms = np.array(n_atoms)

    def _resolve_volume_fractions(self, components, volume_fractions):
        if any([nn <= 0 for nn in volume_fractions]):
            raise ValueError('All volume_fractions must be strictly positive')
        if any([el.density is None for el in components]):
            raise ValueError("All components must have a defined density")
        # Convert volume fractions to mass fractions
        volume_fractions = np.array(volume_fractions)
        volume_fractions /= volume_fractions.sum()
        mass_fractions = np.array([vf * el.density for el, vf in
                                    zip(components, volume_fractions)])
        # Further resolve nested compounds
        self._resolve_mass_fractions(components, mass_fractions)

    def _resolve_molar_fractions(self, components, molar_fractions):
        if any([nn <= 0 for nn in molar_fractions]):
            raise ValueError('All molar_fractions must be strictly positive')
        if any([el.molar_mass is None and el.average_molar_mass is None
                for el in components]):
            raise ValueError("All components must have a defined molar mass")
        # Convert molar fractions to mass fractions
        molar_fractions = np.array(molar_fractions)
        molar_fractions /= molar_fractions.sum()
        molar_masses = np.array([el.molar_mass or el.average_molar_mass
                                 for el in components])
        mass_fractions = np.array([mf * mm for mm, mf in
                                    zip(molar_masses, molar_fractions)])
        self._resolve_mass_fractions(components, mass_fractions)

    def _resolve_mass_fractions(self, components, mass_fractions):
        if any([nn <= 0 for nn in mass_fractions]):
            raise ValueError('All mass_fractions must be strictly positive')
        mass_fractions = np.array(mass_fractions)
        mass_fractions /= mass_fractions.sum()
        while any([el.components is not None for el in components]):
            this_components = []
            this_mass_fractions = []
            for el, fr in zip(components, mass_fractions):
                if el.components is not None:
                    this_components.extend(el.components)
                    this_mass_fractions.extend([ffr * fr for
                                                ffr in el.mass_fractions])
                else:
                    this_components.append(el)
                    this_mass_fractions.append(fr)
            components = this_components
            mass_fractions = this_mass_fractions
        if len(components) != len(mass_fractions):
            raise ValueError('Variables `components` and `mass_fractions` must have the same length')
        # Sum duplicates
        agg = defaultdict(float)
        for el, fr in zip(components, mass_fractions):
            agg[el] += fr
        components, mass_fractions = zip(*agg.items())
        self._components = np.array(components)
        self._mass_fractions = np.array(mass_fractions)
        # Normalise mass fractions
        self._mass_fractions /= self._mass_fractions.sum()

# Cross sections -------------------------------------------------------

    def _get_R(self, A=None):
        """Nuclear radius [m]."""
        A_eff = self.A if A is None else A
        A_eff = np.asarray(A_eff)  # Ensure A_eff is an array
        R = np.where(A_eff > 21,
                 1.1 * A_eff**(1.0/3.0) * 1e-15 * 0.9,
                 1.1 * A_eff**(1.0/3.0) * 1e-15 * 1.05)
        return R

    def _pi_R2_mb(self, A=None):
        """pi*R^2 in mb  (1 fm^2 = 10 mb)."""
        R_fm = self._get_R(A=A) * 1e15
        return np.pi * R_fm**2 * 10.0

    def _compute_glauber_cs(self, A_sig_tot, A_sig_inel, piR2):
        cs_tot   = 2 * piR2 * np.log(1.0 + A_sig_tot / (2 * piR2))
        cs_inel  = piR2 * np.log(1.0 + A_sig_tot / piR2)
        cs_el    = np.maximum(1e-15, cs_tot - cs_inel)

        cs_prod  = piR2 * np.log(1.0 + A_sig_inel / piR2)

        alpha    = A_sig_tot / (2 * piR2 + A_sig_tot)
        cs_sd    = piR2 * (alpha - np.log(1.0 + alpha))

        cs_el_nucleon = np.maximum(1e-15, cs_inel - cs_prod)

        return cs_tot, cs_prod, cs_el, cs_el_nucleon, cs_sd

    def _glauber_element_single(self, sqrt_s, cs_hN, beam, A=None):

        A_eff  = self.A if A is None else A
        piR2   = self._pi_R2_mb(A=A_eff)
        sqrt_s = np.asarray(sqrt_s)

        # cs_hN returns all species, so we need to store all of them.
        A_sig = {
            sp: {
                "tot":  np.zeros_like(sqrt_s),
                "inel": np.zeros_like(sqrt_s),
                "el":   np.zeros_like(sqrt_s),
            }
            for sp in SPECIES
        }

        # Fill hadron–nucleon cross sections
        for i, s in enumerate(sqrt_s):
            values = cs_hN(A_eff, self.Z, s)

            (A_sig["pp"]["tot"][i], A_sig["pp"]["inel"][i], A_sig["pp"]["el"][i],
             A_sig["kmin"]["tot"][i], A_sig["kmin"]["inel"][i], A_sig["kmin"]["el"][i],
             A_sig["kplus"]["tot"][i], A_sig["kplus"]["inel"][i], A_sig["kplus"]["el"][i],
             A_sig["pimin"]["tot"][i], A_sig["pimin"]["inel"][i], A_sig["pimin"]["el"][i],
             A_sig["piplus"]["tot"][i],A_sig["piplus"]["inel"][i],A_sig["piplus"]["el"][i],) = values
        if A_eff < 4:
            return {
                f"cs_{key}_{beam}_GG": A_sig[beam][key]
                for key in ["tot","inel","el"]
            }
        results = {}

        cs_tot, cs_prod, cs_el, cs_el_nucl, cs_sd  = self._compute_glauber_cs(A_sig[beam]["tot"],
                                                                              A_sig[beam]["inel"],
                                                                              piR2)
        results.update({f"cs_tot_{beam}_GG": cs_tot,
                            f"cs_prod_{beam}_GG": cs_prod,
                            f"cs_el_{beam}_GG": cs_el,
                            f"cs_el_nucleon_{beam}_GG": cs_el_nucl,
                            f"cs_sd_{beam}_GG": cs_sd,})
        return results

    def _load_isotopes(self):
        Z = self.Z
        file = Path(__file__).resolve().parent / "isotopes" / f"{Z_TO_SYMBOL.get(Z)}_isotopes.npz"
        data = np.load(file, allow_pickle=True)
        return data

    def _glauber_isotope_to_element(self, log_sqrt_s, cs_hN, species):
        try:
            data = self._load_isotopes()
        except FileNotFoundError:
            return self._glauber_element_single(np.exp(log_sqrt_s), cs_hN, species)

        symbol = Z_TO_SYMBOL.get(self.Z)
        iso_list = ISOTOPES[symbol]["isotopes"]
        stable = [iso for iso in iso_list if iso["abundance"] is not None]
        sp_data = data[species].item()  # the dict for this species
    
        # If self.A matches a specific isotope mass number, use only that one
        matching = [iso for iso in stable if iso["mass_number"] == self.A]
        if matching:
            # Pure isotope — find its index in stable list
            idx = stable.index(matching[0])
            combined = {}
            for k, records in sp_data.items():
                spl = CubicSpline(records[idx]["knots"], records[idx]["y_values"], bc_type='not-a-knot')
                combined[k] = spl(log_sqrt_s)
            return combined

        # Otherwise do abundance-weighted mix as before
        abundance = np.array([iso["abundance"] for iso in stable])
        fraction = abundance / np.sum(abundance)
        combined  = {k: np.zeros_like(log_sqrt_s) for k in sp_data}
        for k, records in sp_data.items():
            for iso_idx, frac in enumerate(fraction):
                spl = CubicSpline(records[iso_idx]["knots"], records[iso_idx]["y_values"], bc_type='not-a-knot')
                combined[k] += frac * spl(log_sqrt_s)
        return combined

    def _glauber_component(self, log_sqrt_s, cs_hN, species):
        """
        Compute component-weighted GG cross sections for a compound material.
        """
        combined = {k: np.zeros_like(log_sqrt_s) 
        for k in self.components[0]._glauber_isotope_to_element(log_sqrt_s, cs_hN, species)}
        for i, comp in enumerate(self.components):
            gg = comp._glauber_isotope_to_element(log_sqrt_s, cs_hN, species)
            for k, v in gg.items():
                combined[k] += self.molar_fractions[i] * v
        return combined

    def _set_attr(self, gg, species, log_sqrt_s, sqrt_s_arr, n_points):
        if species == "pp":
            self._cs_sqrt_s_pp         = sqrt_s_arr
            self._cs_log_sqrt_s_min_pp = float(log_sqrt_s[0])
            self._cs_log_step_pp       = float((log_sqrt_s[-1] - log_sqrt_s[0]) / (n_points - 1))
            self._cs_knots_pp          = log_sqrt_s
            self._n_points_pp          = n_points
            spline_pp = CubicSpline(log_sqrt_s, gg["cs_tot_pp_GG"])
            self._cs_tot_pp_a = spline_pp.c[0]
            self._cs_tot_pp_b = spline_pp.c[1]
            self._cs_tot_pp_c = spline_pp.c[2]
            self._cs_tot_pp_d = spline_pp.c[3]
            del spline_pp
            spline_el_pp = CubicSpline(log_sqrt_s, gg["cs_el_pp_GG"])
            self._cs_el_pp_a = spline_el_pp.c[0]
            self._cs_el_pp_b = spline_el_pp.c[1]
            self._cs_el_pp_c = spline_el_pp.c[2]
            self._cs_el_pp_d = spline_el_pp.c[3]
            del spline_el_pp
            spline_prod_pp = CubicSpline(log_sqrt_s, gg["cs_prod_pp_GG"])
            self._cs_prod_pp_a = spline_prod_pp.c[0]
            self._cs_prod_pp_b = spline_prod_pp.c[1]
            self._cs_prod_pp_c = spline_prod_pp.c[2]
            self._cs_prod_pp_d = spline_prod_pp.c[3]
            del spline_prod_pp
            spline_el_nucleon_pp = CubicSpline(log_sqrt_s, gg["cs_el_nucleon_pp_GG"])
            self._cs_el_nucleon_pp_a = spline_el_nucleon_pp.c[0]
            self._cs_el_nucleon_pp_b = spline_el_nucleon_pp.c[1]
            self._cs_el_nucleon_pp_c = spline_el_nucleon_pp.c[2]
            self._cs_el_nucleon_pp_d = spline_el_nucleon_pp.c[3]
            del spline_el_nucleon_pp
            spline_sd_pp = CubicSpline(log_sqrt_s, gg["cs_sd_pp_GG"])
            self._cs_sd_pp_a = spline_sd_pp.c[0]
            self._cs_sd_pp_b = spline_sd_pp.c[1]
            self._cs_sd_pp_c = spline_sd_pp.c[2]
            self._cs_sd_pp_d = spline_sd_pp.c[3]
            del spline_sd_pp
        elif species == "kmin":
            self._cs_sqrt_s_kmin         = sqrt_s_arr
            self._cs_log_sqrt_s_min_kmin = float(log_sqrt_s[0])
            self._cs_log_step_kmin       = float((log_sqrt_s[-1] - log_sqrt_s[0]) / (n_points - 1))
            self._cs_knots_kmin          = log_sqrt_s
            self._n_points_kmin          = n_points
            spline_kmin = CubicSpline(log_sqrt_s, gg["cs_tot_kmin_GG"])
            self._cs_tot_kmin_a = spline_kmin.c[0]
            self._cs_tot_kmin_b = spline_kmin.c[1]
            self._cs_tot_kmin_c = spline_kmin.c[2]
            self._cs_tot_kmin_d = spline_kmin.c[3]
            del spline_kmin
            spline_el_kmin = CubicSpline(log_sqrt_s, gg["cs_el_kmin_GG"])
            self._cs_el_kmin_a = spline_el_kmin.c[0]
            self._cs_el_kmin_b = spline_el_kmin.c[1]
            self._cs_el_kmin_c = spline_el_kmin.c[2]
            self._cs_el_kmin_d = spline_el_kmin.c[3]
            del spline_el_kmin
            spline_prod_kmin = CubicSpline(log_sqrt_s, gg["cs_prod_kmin_GG"])
            self._cs_prod_kmin_a = spline_prod_kmin.c[0]
            self._cs_prod_kmin_b = spline_prod_kmin.c[1]
            self._cs_prod_kmin_c = spline_prod_kmin.c[2]
            self._cs_prod_kmin_d = spline_prod_kmin.c[3]
            del spline_prod_kmin
            spline_el_nucleon_kmin = CubicSpline(log_sqrt_s, gg["cs_el_nucleon_kmin_GG"])
            self._cs_el_nucleon_kmin_a = spline_el_nucleon_kmin.c[0]
            self._cs_el_nucleon_kmin_b = spline_el_nucleon_kmin.c[1]
            self._cs_el_nucleon_kmin_c = spline_el_nucleon_kmin.c[2]
            self._cs_el_nucleon_kmin_d = spline_el_nucleon_kmin.c[3]
            del spline_el_nucleon_kmin
            spline_sd_kmin = CubicSpline(log_sqrt_s, gg["cs_sd_kmin_GG"])
            self._cs_sd_kmin_a = spline_sd_kmin.c[0]
            self._cs_sd_kmin_b = spline_sd_kmin.c[1]
            self._cs_sd_kmin_c = spline_sd_kmin.c[2]
            self._cs_sd_kmin_d = spline_sd_kmin.c[3]
            del spline_sd_kmin
        elif species == "kplus":
            self._cs_sqrt_s_kplus         = sqrt_s_arr
            self._cs_log_sqrt_s_min_kplus = float(log_sqrt_s[0])
            self._cs_log_step_kplus       = float((log_sqrt_s[-1] - log_sqrt_s[0]) / (n_points - 1))
            self._cs_knots_kplus          = log_sqrt_s
            self._n_points_kplus          = n_points
            spline_kplus = CubicSpline(log_sqrt_s, gg["cs_tot_kplus_GG"])
            self._cs_tot_kplus_a = spline_kplus.c[0]
            self._cs_tot_kplus_b = spline_kplus.c[1]
            self._cs_tot_kplus_c = spline_kplus.c[2]
            self._cs_tot_kplus_d = spline_kplus.c[3]
            del spline_kplus
            spline_el_kplus = CubicSpline(log_sqrt_s, gg["cs_el_kplus_GG"])
            self._cs_el_kplus_a = spline_el_kplus.c[0]
            self._cs_el_kplus_b = spline_el_kplus.c[1]
            self._cs_el_kplus_c = spline_el_kplus.c[2]
            self._cs_el_kplus_d = spline_el_kplus.c[3]
            del spline_el_kplus
            spline_prod_kplus = CubicSpline(log_sqrt_s, gg["cs_prod_kplus_GG"])
            self._cs_prod_kplus_a = spline_prod_kplus.c[0]
            self._cs_prod_kplus_b = spline_prod_kplus.c[1]
            self._cs_prod_kplus_c = spline_prod_kplus.c[2]
            self._cs_prod_kplus_d = spline_prod_kplus.c[3]
            del spline_prod_kplus
            spline_el_nucleon_kplus = CubicSpline(log_sqrt_s, gg["cs_el_nucleon_kplus_GG"])
            self._cs_el_nucleon_kplus_a = spline_el_nucleon_kplus.c[0]
            self._cs_el_nucleon_kplus_b = spline_el_nucleon_kplus.c[1]
            self._cs_el_nucleon_kplus_c = spline_el_nucleon_kplus.c[2]
            self._cs_el_nucleon_kplus_d = spline_el_nucleon_kplus.c[3]
            del spline_el_nucleon_kplus
            spline_sd_kplus = CubicSpline(log_sqrt_s, gg["cs_sd_kplus_GG"])
            self._cs_sd_kplus_a = spline_sd_kplus.c[0]
            self._cs_sd_kplus_b = spline_sd_kplus.c[1]
            self._cs_sd_kplus_c = spline_sd_kplus.c[2]
            self._cs_sd_kplus_d = spline_sd_kplus.c[3]
            del spline_sd_kplus
        elif species == "pimin":
            self._cs_sqrt_s_pimin         = sqrt_s_arr
            self._cs_log_sqrt_s_min_pimin = float(log_sqrt_s[0])
            self._cs_log_step_pimin       = float((log_sqrt_s[-1] - log_sqrt_s[0]) / (n_points - 1))
            self._cs_knots_pimin          = log_sqrt_s
            self._n_points_pimin          = n_points
            spline_pimin = CubicSpline(log_sqrt_s, gg["cs_tot_pimin_GG"])
            self._cs_tot_pimin_a = spline_pimin.c[0]
            self._cs_tot_pimin_b = spline_pimin.c[1]
            self._cs_tot_pimin_c = spline_pimin.c[2]
            self._cs_tot_pimin_d = spline_pimin.c[3]
            del spline_pimin
            spline_el_pimin = CubicSpline(log_sqrt_s, gg["cs_el_pimin_GG"])
            self._cs_el_pimin_a = spline_el_pimin.c[0]
            self._cs_el_pimin_b = spline_el_pimin.c[1]
            self._cs_el_pimin_c = spline_el_pimin.c[2]
            self._cs_el_pimin_d = spline_el_pimin.c[3]
            del spline_el_pimin
            spline_prod_pimin = CubicSpline(log_sqrt_s, gg["cs_prod_pimin_GG"])
            self._cs_prod_pimin_a = spline_prod_pimin.c[0]
            self._cs_prod_pimin_b = spline_prod_pimin.c[1]
            self._cs_prod_pimin_c = spline_prod_pimin.c[2]
            self._cs_prod_pimin_d = spline_prod_pimin.c[3]
            del spline_prod_pimin
            spline_el_nucleon_pimin = CubicSpline(log_sqrt_s, gg["cs_el_nucleon_pimin_GG"])
            self._cs_el_nucleon_pimin_a = spline_el_nucleon_pimin.c[0]
            self._cs_el_nucleon_pimin_b = spline_el_nucleon_pimin.c[1]
            self._cs_el_nucleon_pimin_c = spline_el_nucleon_pimin.c[2]
            self._cs_el_nucleon_pimin_d = spline_el_nucleon_pimin.c[3]
            del spline_el_nucleon_pimin
            spline_sd_pimin = CubicSpline(log_sqrt_s, gg["cs_sd_pimin_GG"])
            self._cs_sd_pimin_a = spline_sd_pimin.c[0]
            self._cs_sd_pimin_b = spline_sd_pimin.c[1]
            self._cs_sd_pimin_c = spline_sd_pimin.c[2]
            self._cs_sd_pimin_d = spline_sd_pimin.c[3]
            del spline_sd_pimin
        elif species == "piplus":
            self._cs_sqrt_s_piplus         = sqrt_s_arr
            self._cs_log_sqrt_s_min_piplus = float(log_sqrt_s[0])
            self._cs_log_step_piplus       = float((log_sqrt_s[-1] - log_sqrt_s[0]) / (n_points - 1))
            self._cs_knots_piplus          = log_sqrt_s
            self._n_points_piplus          = n_points
            spline_piplus = CubicSpline(log_sqrt_s, gg["cs_tot_piplus_GG"])
            self._cs_tot_piplus_a = spline_piplus.c[0]
            self._cs_tot_piplus_b = spline_piplus.c[1]
            self._cs_tot_piplus_c = spline_piplus.c[2]
            self._cs_tot_piplus_d = spline_piplus.c[3]
            del spline_piplus
            spline_el_piplus = CubicSpline(log_sqrt_s, gg["cs_el_piplus_GG"])
            self._cs_el_piplus_a = spline_el_piplus.c[0]
            self._cs_el_piplus_b = spline_el_piplus.c[1]
            self._cs_el_piplus_c = spline_el_piplus.c[2]
            self._cs_el_piplus_d = spline_el_piplus.c[3]
            del spline_el_piplus
            spline_prod_piplus = CubicSpline(log_sqrt_s, gg["cs_prod_piplus_GG"])
            self._cs_prod_piplus_a = spline_prod_piplus.c[0]
            self._cs_prod_piplus_b = spline_prod_piplus.c[1]
            self._cs_prod_piplus_c = spline_prod_piplus.c[2]
            self._cs_prod_piplus_d = spline_prod_piplus.c[3]
            del spline_prod_piplus
            spline_el_nucleon_piplus = CubicSpline(log_sqrt_s, gg["cs_el_nucleon_piplus_GG"])
            self._cs_el_nucleon_piplus_a = spline_el_nucleon_piplus.c[0]
            self._cs_el_nucleon_piplus_b = spline_el_nucleon_piplus.c[1]
            self._cs_el_nucleon_piplus_c = spline_el_nucleon_piplus.c[2]
            self._cs_el_nucleon_piplus_d = spline_el_nucleon_piplus.c[3]
            del spline_el_nucleon_piplus
            spline_sd_piplus = CubicSpline(log_sqrt_s, gg["cs_sd_piplus_GG"])
            self._cs_sd_piplus_a = spline_sd_piplus.c[0]
            self._cs_sd_piplus_b = spline_sd_piplus.c[1]
            self._cs_sd_piplus_c = spline_sd_piplus.c[2]
            self._cs_sd_piplus_d = spline_sd_piplus.c[3]
            del spline_sd_piplus
        else:
            raise ValueError(f"Unknown species {species}")

    def _get_material_cross_sections(self, n_points=None, sqrt_s_min=None, sqrt_s_max=None):
        splines, grid_min, grid_max = load_all_splines()
        _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_, cs_hN = make_nucleon_cs(splines)
        n_points = n_points or N_CS_POINTS

        if self._components is None:
            if self.A <= 62:
                self._nuclear_slope      = 14.5 * self.A**(2./3.)
                self._nuclear_slope_pion = 14.5 * self.A**(2./3.)
            else:
                self._nuclear_slope      = 60.0 * self.A**(1./3.)
                self._nuclear_slope_pion = 60.0 * self.A**(1./3.)
        else:
            nuclear_slope_14 = sum(self.molar_fractions[i] * comp.A**(2./3.)
                                   for i, comp in enumerate(self.components) if comp.A <= 62) * 14.5
            nuclear_slope_60 = sum(self.molar_fractions[i] * comp.A**(1./3.)
                                   for i, comp in enumerate(self.components) if comp.A > 62) * 60.0
            self._nuclear_slope      = nuclear_slope_14 + nuclear_slope_60
            self._nuclear_slope_pion = nuclear_slope_14 + nuclear_slope_60
        for beam in SPECIES:
            smin = sqrt_s_min or grid_min[beam]
            smax = sqrt_s_max or grid_max[beam]
            sqrt_s_arr = np.logspace(np.log10(smin), np.log10(smax), n_points)
            log_sqrt_s = np.log(sqrt_s_arr)

            if self._components is None:
                gg = self._glauber_isotope_to_element(log_sqrt_s, cs_hN, beam)
            else:
                gg = self._glauber_component(log_sqrt_s, cs_hN, beam)
            self._set_attr(gg, beam, log_sqrt_s, sqrt_s_arr, n_points)

    def evaluate_glauber_spline(self, sqrt_s, key, particle_id):
        # Key: 0 = total, 1 = inelastic, 2 = elastic, 3 = production, 4 = single diffractive
        self.compile_kernels(only_if_needed=True)
        return self._context.kernels.evaluate_glauber_spline(
            material=self,
            sqrt_s=float(sqrt_s),
            key=int(key),
            particle_id=int(particle_id)
        )
    # ===========
    # === API ===
    # ===========

    @classmethod
    def from_dict(cls, dct):
        from xcoll.materials.database import db as mdb
        this_cls = dct.get('__class__', cls.__name__)
        if cls.__name__ != this_cls:
            if this_cls == 'Material':
                return Material.from_dict(dct)
            elif this_cls == 'RefMaterial':
                return RefMaterial.from_dict(dct)
            else:
                raise ValueError(f"Unknown material class {this_cls} in from_dict")
        # Create instance but do not set names (to avoid syncing with the database)
        name = dct.pop('name', None)
        short_name = dct.pop('short_name', None)
        fluka_name = dct.pop('fluka_name', None)
        geant4_name = dct.pop('geant4_name', None)
        mat = super().from_dict(dct)
        if name in mdb:
            # Set names without storing in database
            mat._name = name
            mat._short_name = short_name
            mat._fluka_name = fluka_name
            mat._geant4_name = geant4_name
            if mdb[name] == mat:
                # Return the one from the database
                return mdb[name]
            else:
                # Return the new one but mark as out-of-sync
                print(f"Warning: Material {name} is out of sync with database.")
                mat._out_of_sync = True
                return mat
        else:
            # Store in database
            mat.name = name
            mat.short_name = short_name
            mat.fluka_name = fluka_name
            mat.geant4_name = geant4_name
            return mat

    def to_dict(self, *args, **kwargs):
        dct = super().to_dict(*args, **kwargs)
        if not self._radiation_length_set_manually:
            dct.pop('radiation_length', None)
        if not self._excitation_energy_set_manually:
            dct.pop('excitation_energy', None)
        if self.n_atoms is not None:
            dct.pop('mass_fractions', None)
        if 'components' in dct:
            dct['components'] = [el.to_dict() for el in self.components]
        return dct

    def adapt(self, inplace=False, **kwargs):
        thisdict = self.to_dict()
        thisdict.pop('__class__')
        thisdict.pop('name', None)
        thisdict.pop('short_name', None)
        thisdict.pop('fluka_name', None)
        thisdict.pop('geant4_name', None)
        thisdict.pop('info', None)
        thisdict.update(kwargs)
        if inplace:
            for kk, vv in thisdict.items():
                if kk in ['A', 'Z', 'components', 'n_atoms', 'mass_fractions',
                          'volume_fractions', 'molar_fractions',
                          'atomic_fractions']:
                    if kk in kwargs:
                        raise ValueError(f"Cannot adapt {kk} inplace")
                    continue
                setattr(self, kk, vv)
            return self
        return self.__class__(**thisdict)

    def invalidate(self):
        self.Z = None
        self.A = None
        self._components = None
        self._n_atoms = None
        self._mass_fractions = None
        self.density = None
        self.name = None
        self.short_name = None
        self.fluka_name = None
        self.geant4_name = None

    def __hash__(self):
        if self.name is None:
            raise TypeError('Material must have a name to be hashable.')
        return hash(self.name)

    def __repr__(self):
        return f"<{self} at {hex(id(self))}>"

    def __str__(self):
        cls = self.__class__.__name__
        if self.is_elemental:
            typ = 'Elemental'
            Z = int(self.Z)
            if not np.isclose(self.A, int(self.A), atol=1e-15):
                A = f'{self.A:.3f}'
            else:
                A = f'[{int(self.A)}]'
            comp = f"Z={Z}, A={A}"
        elif self.is_compound:
            typ = 'Compound'
            comp = [f"{nn}: {round(fr, 2):g}" for nn, fr in self.composition]
            comp = f"[{', '.join(comp)}]"
        elif self.is_mixture:
            typ = 'Mixture'
            comp = [f"{nn}: {round(fr, 5):g}" for nn, fr in self.composition]
            comp = f"[{', '.join(comp)}]"
        else:
            return f"Invalid {cls}({self.name})"
        comp += f", density={self.density:.4f} g/cm^3"
        if self.full_everest_crystal_supported:
            everest = ' [full support for Everest crystals]'
        elif self.full_everest_supported:
            everest = ' [full support for Everest]'
        else:
            everest = ''
        return f"{typ} {cls}({self.name}, {comp}){everest}"

    def __eq__(self, other):
        if isinstance(other, str):
            from xcoll.materials.database import db as mdb
            try:
                other = mdb[other]
            except KeyError:
                return False
        if not isinstance(other, Material):
            return False
        keys = set(self.to_dict().keys())
        if keys != set(other.to_dict().keys()):
            return False
        dct_self = self.to_dict().copy()
        dct_other = other.to_dict().copy()
        if 'components' in dct_self:
            if 'components' not in dct_other:
                return False
            if len(dct_self['components']) != len(dct_other['components']):
                return False
            for el1, el2 in zip(dct_self['components'], dct_other['components']):
                if not deep_equal(el1, el2):
                    return False
            dct_self.pop('components')
            dct_other.pop('components')
        return deep_equal(dct_self, dct_other)

    def __setattr__(self, name, value):
        if name not in ['_xobject', '_frozen']:
            self._assert_not_frozen(name)
        return super().__setattr__(name, value)

    def _assert_not_frozen(self, name):
        if hasattr(self, '_frozen') and self._frozen:
            raise ValueError(f'Material is frozen, cannot change attribute {name}.')

    @property
    def full_everest_supported(self):
        if self.excitation_energy is None:
            return False
        if self.nuclear_radius is None:
            return False
        if self.nuclear_elastic_slope is None:
            return False
        if self.hcut is None:
            return False
        if self.cross_section is None:
            return False
        return True

    @property
    def full_everest_crystal_supported(self):
        if not self.full_everest_supported:
            return False
        if self.crystal_plane_distance is None:
            return False
        if self.crystal_potential is None:
            return False
        if self.nuclear_collision_length is None:
            return False
        if self.eta is None:
            return False
        return True

    @property
    def is_elemental(self):
        return self.A is not None

    @property
    def is_compound(self):
        return self.components is not None and self.n_atoms is not None

    @property
    def is_mixture(self):
        return self.components is not None and self.n_atoms is None

    @property
    def composition(self):
        if self.n_atoms is not None:
            return [[el.short_name or el.name, nn.tolist()]
                    for el, nn in zip(self.components, self.n_atoms)]
        elif self.mass_fractions is not None:
            return [[el.short_name or el.name, fr.tolist()]
                    for el, fr in zip(self.components, self.mass_fractions)]

    def update_vars(self):
        # Update average Z over A
        if self.components is not None:
            self._ZA_mean = _average_Z_over_A(self.components,
                                              self.mass_fractions)
        elif self.Z is not None and self.A is not None:
            self._ZA_mean = self.Z / self.A
        else:
            self._ZA_mean = -1
        # Update effective Z
        if self.components is not None:
            self._Z2_eff = _effective_Z2(self.components,
                                         self.molar_fractions)
        elif self.Z is not None:
            self._Z2_eff = self.Z*self.Z
        else:
            self._Z2_eff = -1
        # Update _atoms_per_volume
        if self.atoms_per_volume is not None:
            self._atoms_per_volume = self.atoms_per_volume
        else:
            self._atoms_per_volume = -1
        # Update _num_nucleons_eff
        if self.num_effective_nucleons is not None:
            self._num_nucleons_eff = self.num_effective_nucleons
        else:
            self._num_nucleons_eff = -1
        # Update radiation length and excitation energy if not set manually
        if not self._radiation_length_set_manually:
            self.radiation_length = None   # Recompute radiation length
        if not self._excitation_energy_set_manually:
            self.excitation_energy = None  # Recompute excitation energy

    # =======================
    # === Main Properties ===
    # =======================

    @property
    def Z(self):
        return self._Z

    @Z.setter
    def Z(self, val):
        self._assert_not_frozen('Z')
        if self.is_compound:
            raise ValueError('Cannot set Z for a compound material')
        elif self.is_mixture:
            raise ValueError('Cannot set Z for a mixture material')
        if val is None:
            self._Z = None
        else:
            if val <= 0:
                raise ValueError('Z must be strictly positive')
            self._Z = int(round(val))
        self.update_vars()

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, val):
        self._assert_not_frozen('A')
        if self.is_compound:
            raise ValueError('Cannot set A for a compound material')
        elif self.is_mixture:
            raise ValueError('Cannot set A for a mixture material')
        if val is None:
            self._A = None
        else:
            if val <= 0:
                raise ValueError('A must be strictly positive')
            self._A = val
        self.update_vars()

    @property
    def components(self):
        if self._components is not None:
            return np.array(self._components)

    @property
    def n_atoms(self):
        if self._n_atoms is not None:
            return np.array(self._n_atoms)

    @property
    def mass_fractions(self):
        if self.n_atoms is not None:
            fr = np.array([nn * el.A for el, nn
                                in zip(self.components, self.n_atoms)])
            return fr / fr.sum()
        elif self._mass_fractions is not None:
            return np.array(self._mass_fractions)

    @property
    def molar_fractions(self):
        if self.n_atoms is not None:
            return self.n_atoms / sum(self.n_atoms)
        elif self.components is not None:
            # Use molar_mass such that this can be used for nested compounds
            mf = np.array([fr / el.molar_mass for el, fr
                                in zip(self.components, self.mass_fractions)])
            return mf / mf.sum()

    @property
    def atomic_fractions(self):
        if self.components is not None:
            return self.molar_fractions

    @property
    def volume_fractions(self):
        if self.components is not None:
            vf = np.array([fr / el.density for el, fr
                                    in zip(self.components, self.mass_fractions)])
            return vf / vf.sum()

    @property
    def molar_mass(self):
        if self.n_atoms is not None:
            return sum([nn * el.A for el, nn
                                in zip(self.components, self.n_atoms)])
        else:
            return self.A

    @property
    def average_molar_mass(self):
        if self._mass_fractions is not None:
            # Use molar_mass such that this can be used for nested compounds
            return sum([mf * el.molar_mass for el, mf
                               in zip(self.components, self.molar_fractions)])

    @property
    def density(self):
        if self._density > 0:
            return self._density

    @density.setter
    def density(self, val):
        self._assert_not_frozen('density')
        if val is None:
            self._density = -1
        else:
            if val <= 0:
                raise ValueError('density must be strictly positive')
            self._density = val
        self.update_vars()

    @property
    def radiation_length(self):
        if self._radiation_length > 0:
            return self._radiation_length

    @radiation_length.setter
    def radiation_length(self, val):
        self._assert_not_frozen('radiation_length')
        if val is None:
            if self.components is not None and self.mass_fractions is not None \
            and self.density is not None:
                self._radiation_length = _combine_radiation_lengths(
                            self.components, self.mass_fractions, self.density)
            elif self.Z is not None and self.A is not None \
            and self.density is not None:
                self._radiation_length = _approximate_radiation_length(
                                                self.Z, self.A, self.density)
            else:
                self._radiation_length = -1
            self._radiation_length_set_manually = False
        else:
            if val <= 0:
                raise ValueError('radiation_length must be strictly positive')
            self._radiation_length = val
            self._radiation_length_set_manually = True

    @property
    def excitation_energy(self):
        if self._excitation_energy > 0:
            return self._excitation_energy

    @excitation_energy.setter
    def excitation_energy(self, val):
        self._assert_not_frozen('excitation_energy')
        if val is None:
            if self.components is not None:
                self._excitation_energy = _combine_excitation_energies(
                                        self.components, self.mass_fractions)
            else:
                self._excitation_energy = _default_excitation_energies.get(
                                                                self.Z, -1.)
            self._excitation_energy_set_manually = False
        else:
            if val <= 0:
                raise ValueError('excitation_energy must be strictly positive')
            self._excitation_energy = val
            self._excitation_energy_set_manually = True

    @property
    def nuclear_radius(self):
        if self._nuclear_radius > 0:
            return self._nuclear_radius

    @nuclear_radius.setter
    def nuclear_radius(self, val):
        self._assert_not_frozen('nuclear_radius')
        if val is None:
            self._nuclear_radius = -1
        else:
            if val <= 0:
                raise ValueError('nuclear_radius must be strictly positive')
            self._nuclear_radius = val

    @property
    def nuclear_elastic_slope(self):
        if self._nuclear_elastic_slope > 0:
            return self._nuclear_elastic_slope

    @nuclear_elastic_slope.setter
    def nuclear_elastic_slope(self, val):
        self._assert_not_frozen('nuclear_elastic_slope')
        if val is None:
            self._nuclear_elastic_slope = -1
        else:
            if val <= 0:
                raise ValueError('nuclear_elastic_slope must be strictly positive')
            self._nuclear_elastic_slope = val

    @property
    def cross_section(self):
        if self._cross_section[0] >= 0:
            return self._cross_section

    @cross_section.setter
    def cross_section(self, val):
        self._assert_not_frozen('cross_section')
        if val is None:
            self._cross_section = [-1., -1., -1., -1., -1., -1.]
        else:
            if not hasattr(val, '__iter__') or isinstance(val, str) or len(val) != 6:
                raise ValueError('cross_section must be a list of 6 values')
            for v in val:
                if v < 0:
                    raise ValueError('cross_section values must be positive')
            self._cross_section = val

    @property
    def hcut(self):
        if self._hcut > 0:
            return self._hcut

    @hcut.setter
    def hcut(self, val):
        self._assert_not_frozen('hcut')
        if val is None:
            self._hcut = -1
        else:
            if val <= 0:
                raise ValueError('hcut must be strictly positive')
            self._hcut = val

    @property
    def crystal_plane_distance(self):
        if self._crystal_plane_distance > 0:
            return self._crystal_plane_distance

    @crystal_plane_distance.setter
    def crystal_plane_distance(self, val):
        self._assert_not_frozen('crystal_plane_distance')
        if val is None:
            self._crystal_plane_distance = -1
        elif val <= 0:
            raise ValueError('`crystal_plane_distance` must be strictly positive')
        self._crystal_plane_distance = val

    @property
    def crystal_potential(self):
        if self._crystal_potential > 0:
            return self._crystal_potential

    @crystal_potential.setter
    def crystal_potential(self, val):
        self._assert_not_frozen('crystal_potential')
        if val is None:
            self._crystal_potential = -1
        elif val <= 0:
            raise ValueError('`crystal_potential` must be strictly positive')
        self._crystal_potential = val

    @property
    def nuclear_collision_length(self):
        if self._nuclear_collision_length > 0:
            return self._nuclear_collision_length

    @nuclear_collision_length.setter
    def nuclear_collision_length(self, val):
        self._assert_not_frozen('nuclear_collision_length')
        if val is None:
            self._nuclear_collision_length = -1
        elif val <= 0:
            raise ValueError('`nuclear_collision_length` must be strictly positive')
        self._nuclear_collision_length = val

    @property
    def eta(self):
        if self._eta > 0:
            return self._eta

    @eta.setter
    def eta(self, val):
        self._assert_not_frozen('eta')
        if val is None:
            self._eta = -1
        elif val <= 0:
            raise ValueError('`eta` must be strictly positive')
        self._eta = val

    @property
    def electron_density(self):
        if self._ZA_mean > 0:
            return self._ZA_mean * sc.Avogadro # [electrons/g]

    @property
    def plasma_energy(self):
        if self._ZA_mean > 0 and self.density is not None:
            return np.sqrt(self._ZA_mean * self.density) * 28.816

    @property
    def atoms_per_volume(self):
        # TODO: do we want to use a different A? Maybe an average instead of molar mass?
        mass = self.molar_mass or self.average_molar_mass
        if self.density is not None and mass is not None:
            return self.density * 1.e6 * sc.Avogadro / mass # [atoms/m^3]

    @property
    def num_effective_nucleons(self):
        # See thesis Nuria (3.14)
        dn = 0.415 # Width of outer ring of nucleus that participates in nuclear interactions
        # TODO: do we want to use a different A? Maybe an average instead of molar mass?
        A1_3 = None
        if self.A is not None:
            A1_3 = self.A**(1./3.)
        elif self.components is not None:
            A1_3 = sum([mf * (el.A**(1./3.)) for el, mf
                            in zip(self.components, self.molar_fractions)])
        if A1_3 is not None:
            return 2*np.pi*dn*(3/4/np.pi)**(1/3) * A1_3

    @property
    def isotopes(self):
        return self._isotopes

    @isotopes.setter
    def isotopes(self, val):
        if val is None:
            self._isotopes = None
            return

        if isinstance(val, str):
            if val not in ISOTOPES:
                raise ValueError(f"Unknown element symbol '{val}' for isotopes.")
            val = ISOTOPES[val].get('isotopes', [])
        elif isinstance(val, dict):
            val = val.get('isotopes', [val])

        if not hasattr(val, '__iter__') or isinstance(val, (str, bytes)):
            raise ValueError("`isotopes` must be an iterable, dict, symbol string, or None.")

        isotopes = []
        for iso in val:
            if isinstance(iso, IsotopeData):
                this_iso = iso
            elif isinstance(iso, dict):
                this_iso = IsotopeData(
                    mass_number=int(iso['mass_number']),
                    atomic_mass=float(iso['atomic_mass']),
                    abundance=None if iso.get('abundance', None) is None else float(iso['abundance']),
                )
            else:
                raise ValueError("Each isotope must be IsotopeData or a dict.")

            if this_iso.abundance is not None:
                isotopes.append(this_iso)

        self._isotopes = tuple(isotopes) if isotopes else None

    # =======================
    # === Meta Properties ===
    # =======================

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, val):
        self._assert_not_frozen('state')
        if val is None:
            self._state = None
        elif val.lower() not in ['solid', 'liquid', 'gas']:
            raise ValueError("Attribute `state` must be one of 'solid', "
                             "'liquid', 'gas', or None")
        else:
            self._state = val.lower()

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, val):
        self._assert_not_frozen('temperature')
        if val is None:
            self._temperature = None
        elif val <= 0:
            raise ValueError('Temperature must be in Kelvin and strictly positive')
        else:
            self._temperature = val

    @property
    def pressure(self):
        # In atm = 101.325 kPa
        return self._pressure

    @pressure.setter
    def pressure(self, val):
        self._assert_not_frozen('pressure')
        if val is None:
            self._pressure = None
        elif val <= 0:
            raise ValueError('Pressure must be in Pascal and strictly positive')
        else:
            self._pressure = val

    @property
    def info(self):
        return self._info

    @info.setter
    def info(self, info):
        self._assert_not_frozen('info')
        self._info = info


    # =======================
    # === Name Properties ===
    # =======================

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._assert_not_frozen('name')
        old_name = self.name
        self._name = name.strip() if name is not None else None
        # Sync database
        from xcoll.materials.database import db as mdb
        if name is None:
            if old_name is not None:
                mdb.remove_material(old_name)
        else:
            if old_name is None:
                mdb[name] = self
            elif old_name == name:
                mdb.update_material(name, self)
            else:
                mdb.rename_material(old_name, name)

    @property
    def short_name(self):
        return self._short_name

    @short_name.setter
    def short_name(self, name):
        self._assert_not_frozen('short_name')
        if name is not None and self.name is None:
            raise ValueError('Material must have a name to set `short_name`')
        old_name = self.short_name
        self._short_name = name.strip() if name is not None else None
        # Sync database
        from xcoll.materials.database import db as mdb
        if name is None:
            if old_name is not None:
                mdb.remove_alias(old_name)
        else:
            if old_name is None:
                mdb[name] = self
            else:
                mdb.rename_alias(old_name, name)

    @property
    def fluka_name(self):
        return self._fluka_name

    @fluka_name.setter
    def fluka_name(self, name):
        self._assert_not_frozen('fluka_name')
        old_name = self.fluka_name
        self._fluka_name = name.strip()[:8] if name is not None else None
        # Sync database
        from xcoll.materials.database import db as mdb
        if name is None:
            if old_name is not None:
                mdb.fluka.remove_alias(old_name)
        else:
            if old_name is None:
                mdb.fluka[name] = self
            else:
                mdb.fluka.rename_alias(old_name, name)

    @property
    def geant4_name(self):
        return self._geant4_name

    @geant4_name.setter
    def geant4_name(self, name):
        self._assert_not_frozen('geant4_name')
        old_name = self.geant4_name
        self._geant4_name = name.strip() if name is not None else None
        # Sync database
        from xcoll.materials.database import db as mdb
        if name is None:
            if old_name is not None:
                mdb.geant4.remove_alias(old_name)
        else:
            if old_name is None:
                mdb.geant4[name] = self
            else:
                mdb.geant4.rename_alias(old_name, name)


    # =======================
    # === Code Generation ===
    # =======================

# The largest atomic number that can be handled by FLUKA is 100.
# *...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+...
# MATERIAL        AAAA                BBBB                          CCCCHYDROGEN
# MAT-PROP        DDDD                EEEE  HYDROGEN
# AAAA: Atomic number (0 for compounds/mixtures
# BBBB: density in g/cm3
# CCCCL: optional, specify isotope (A) but only integers allowed
# DDDD: optional, pressure if gas in atm
# EEEE: optional, ionisation potential in eV
#
# Composition of compounds and mixtures is defined in FLUKA input files:
# Negative values -> mass fractions
# Positive values -> number of atoms / atomic ratios
# * Antico
# * ..+....1....+....2....+....3....+....4....+....5....+....6....+....7....+..
# MATERIAL                             2.7                              ANTICO
# COMPOUND      -97.25  ALUMINUM       -.6   SILICON      -0.5      IRONANTICO
# COMPOUND         -.4  MAGNESIU       -.4  MANGANES      -.35  CHROMIUMANTICO
# COMPOUND         -.2  TITANIUM       -.2      ZINC       -.1    COPPERANTICO

    def _generate_fluka_code(self):
        if self.fluka_name and self.fluka_name.startswith('XCOLL') \
        and self._generated_fluka_code is not None:
            return # Already generated
        if self.fluka_name and not self.fluka_name.startswith('XCOLL'):
            raise ValueError(f'Material already known to FLUKA ({self.fluka_name}).')
        if self.name is None:
            raise ValueError('Material must have a name to generate fluka code.')
        if self.components is not None:
            for mat in self.components:
                if mat.fluka_name is None:
                    mat._generate_fluka_code()
        # Create unique name
        from xcoll.materials.database import db as mdb
        existing_ids = [int(nn[5:]) for nn in mdb.fluka.keys() if nn.startswith('XCOLL')]
        max_id = max(existing_ids) if existing_ids else 0
        new_id = min(set(range(max_id+2)) - set(existing_ids))
        name = f'XCOLL{new_id:03d}'
        # Define single element material
        from xcoll.scattering_routines.fluka.environment import format_fluka_float
        if self.components is None:
            if self.Z > 100:
                raise ValueError('FLUKA cannot handle elements with Z > 100.')
            P = format_fluka_float(self.pressure) if self.pressure else 10*' '
            exc = format_fluka_float(self.excitation_energy) if self.excitation_energy else 10*' '
            code = f"""\
*...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+...
MATERIAL  {format_fluka_float(self.Z)}          {format_fluka_float(self.density)}                              {name}
MAT-PROP  {P}          {exc}  {name}
"""
        else:
            code = f"""\
*...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+...
MATERIAL                      {format_fluka_float(self.density)}                              {name}
"""
            if self.n_atoms is not None:
                frac = [self.n_atoms[i:i + 3] for i in range(0, len(self.n_atoms), 3)]
                comp = [self.components[i:i + 3] for i in range(0, len(self.components), 3)]
                for ff, cc in zip(frac, comp):
                    line = "COMPOUND  "
                    for i, (fff, ccc) in enumerate(zip(ff, cc)):
                        line += format_fluka_float(fff)
                        line += f"{ccc.fluka_name[:10]:>10}"
                    line += (2-i)*20*' '
                    line += f"{name}"
                    code += line + "\n"
            else:
                frac = [self.mass_fractions[i:i + 3] for i in range(0, len(self.mass_fractions), 3)]
                comp = [self.components[i:i + 3] for i in range(0, len(self.components), 3)]
                for ff, cc in zip(frac, comp):
                    line = "COMPOUND  "
                    for i, (fff, ccc) in enumerate(zip(ff, cc)):
                        line += format_fluka_float(-fff)
                        line += f"{ccc.fluka_name[:10]:>10}"
                    line += (2-i)*20*' '
                    line += f"{name}"
                    code += line + "\n"
        frozen = self._frozen
        self._frozen = False
        self.fluka_name = name
        self._generated_fluka_code = code
        self._frozen = frozen

    def _generate_geant4_code(self):
        if self.geant4_name and self.geant4_name.startswith('Xcoll') \
        and self._generated_geant4_code is not None:
            return # Already generated
        if self.geant4_name and not self.geant4_name.startswith('Xcoll'):
            raise ValueError(f'Material already has a Geant4 name assigned: {self.geant4_name}.')
        if self.name is None:
            raise ValueError('Material must have a name to generate Geant4 code.')
        if self.components is not None:
            if any([el.geant4_name is None for el in self.components]):
                raise ValueError('All components must have a geant4_name to generate Geant4 code.')
        name = f'Xcoll{self.name}'
        if self.components is not None:
            code = f"{name} : matdef, density={self.density}"
        else:
            code = f"{name} : matdef, Z={self.Z}, A={self.A}, density={self.density}"
        if self.temperature:
            code += f", T={self.temperature}"
        if self.pressure:
            code += f", P={self.pressure}"
        if self.state:
            code += f', state="{self.state}"'
        if self.components is not None:
            components = [f'"{el.geant4_name}"' for el in self.components]
            code += f", components=[{','.join(components)}]"
            if self.n_atoms is not None and np.all([int(nn) == nn for nn in self.n_atoms]):
                # Geant4 wants integer numbers of atoms
                code += f", componentsWeights={{{','.join([f'{nn}' for nn in self.n_atoms])}}};"
            else:
                code += f", componentsFractions={{{','.join([f'{nn}' for nn in self.mass_fractions])}}}"
        code += ";"
        frozen = self._frozen
        self._frozen = False
        self.geant4_name = name
        self._generated_geant4_code = code
        self._frozen = frozen


class RefMaterial(Material):
    ''' A string to pass a material that exists in Geant4 or FLUKA.'''
    _xofields = Material._xofields
    _depends_on = [Material]
    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            if 'fluka_name' not in kwargs and 'geant4_name' not in kwargs:
                raise ValueError("ReferenceMaterial must have a fluka_name or geant4_name.")
            fluka_name = kwargs.pop('fluka_name', None)
            geant4_name = kwargs.pop('geant4_name', None)
            if fluka_name is not None and geant4_name is not None:
                raise ValueError("ReferenceMaterial can have only one of fluka_name or geant4_name.")
            # Fake values to pass initialisation
            kwargs['Z'] = 1
            kwargs['A'] = 1
            kwargs['density'] = 1
        super().__init__(**kwargs)
        # Invalidate all properties
        self.invalidate()
        # Set name only now to avoid syncing with database
        self._name = fluka_name or geant4_name
        if fluka_name is not None:
            self._fluka_name = fluka_name
        if geant4_name is not None:
            self._geant4_name = geant4_name
        self._frozen = True  # ReferenceMaterial is always frozen

    def __str__(self):
        return f"RefMaterial '{self.name}'"


_DEFAULT_MATERIAL = Material(Z=1, A=1, density=1)
_DEFAULT_MATERIAL.invalidate()


def _resolve_material(material, allow_none=None, ref=None, everest_crystal=False):
    if material is None:
        if allow_none is None:
            return _DEFAULT_MATERIAL
        elif allow_none:
            return None
        raise ValueError('Material cannot be None!')
    elif isinstance(material, dict):
        material = Material.from_dict(material)
    elif isinstance(material, str):
        from xcoll.materials.database import db as mdb
        material = mdb[material]
    elif ref is not None and isinstance(material, RefMaterial):
        if getattr(material, f'{ref}_name', None) is None:
            raise ValueError(f"RefMaterial {material} does not have a {ref} name!")
    if not isinstance(material, Material):
        raise ValueError(f"Invalid material of type {type(material)}!")
    if everest_crystal and not material.full_everest_crystal_supported:
        raise ValueError(f"Material {material.name} does not have full Everest crystal support!")
    return material
