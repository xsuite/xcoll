# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from unicodedata import name

import numpy as np
from numbers import Number
from contextlib import contextmanager

import xtrack as xt
import xtrack.particles.pdg as pdg

from ..pretty_print import style


class PhysicsSettingsHelper:
    """Helper class to manage physics settings for scattering routines.
    """

    # Prepend 'return_' to get the property name for the return flags
    _return_flags = {'photons': [],
                     'leptons': ['electrons', 'muons', 'tauons', 'neutrinos'],
                     'baryons': ['protons', 'neutrons', 'other_baryons'],
                     'mesons':  ['pions', 'kaons', 'other_mesons'],
                     'ions':    []}
    _global_return_flags = ['all', 'all_charged', 'none']
    _return_modifiers = ['neutral']

    # Append '_cut' to get the property name for cut definitions
    _cut_definitions = ['hadron_lower_momentum', 'photon_lower_momentum',
                        'electron_lower_momentum', 'relative_energy']

    # Prepend 'include_' to get the property name for the physics processes
    _include_processes = ['showers', 'single_coulomb', 'multiple_coulomb',
                          'ionisation_fluctuations', 'pair_production',
                          'bremsstrahlung', 'elastic', 'inelastic']


    # ===========
    # === API ===
    # ===========

    def __init__(self, engine):
        from xcoll.scattering_routines.engine import BaseEngine
        if not isinstance(engine, BaseEngine):
            raise ValueError("`engine` has to be an instance of BaseEngine!")
        with self.__class__._in_constructor(self):
            self._engine = engine
            self.reset()

    def show(self, full=False):
        """Print the physics settings."""
        print(self._str(format=True))

    @property
    def all_flags(self):
        all_flags  = self._global_return_flags.copy()
        all_flags += self._return_flags.keys()
        all_flags += self._return_modifiers
        all_flags += [fff for ff in self._return_flags.values() for fff in ff]
        result  = [f'return_{ff}' for ff in all_flags]
        result += [f'{ff}_cut' for ff in self._cut_definitions]
        result += [f'include_{pp}' for pp in self._include_processes]
        return result

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))} (use .show() " \
             + f"to see the contents)>"

    def __str__(self):
        return self._str(format=False)

    def _print(self, *args, **kwargs):
        return self._engine._print(*args, **kwargs)

    def _error(self, mess, error_class=ValueError):
        self._engine.stop()
        raise error_class(mess)

    def _str(self, format):
        final_message = ''
        veto = self._engine._physics_settings_veto_list

        # Global return flags
        mess = ''
        flags = [f'return_{ff}' for ff in self._global_return_flags
                 if f'return_{ff}' not in veto]
        for i, flag in enumerate(flags):
            prefix = "└" if i == len(flags) - 1 else "├"
            val = getattr(self, flag)
            name = f"{flag.replace('return_', '').replace('_', ' ')}:"
            mess += f"  {prefix} {name:15} {val}\n"
        if mess != '':
            title = "Global return flags:"
            title = style(f"{title:25}", bold=True, colour='forest_green',
                          enabled=format)
            title += style(" (prepend 'return_' to get the property name)",
                           dim=True, italic=True, colour='forest_green',
                           enabled=format)
            final_message += f"{title}\n{mess}\n"

        # Global return modifiers
        mess = ''
        flags = [f'return_{ff}' for ff in self._return_modifiers
                 if f'return_{ff}' not in veto]
        for i, flag in enumerate(flags):
            prefix = "└" if i == len(flags) - 1 else "├"
            val = getattr(self, flag)
            name = f"{flag.replace('return_', '').replace('_', ' ')}:"
            mess += f"  {prefix} {name:15} {val}\n"
        if mess != '':
            title = "Global return modifiers:"
            title = style(f"{title:25}", bold=True, colour='forest_green',
                          enabled=format)
            title += style(" (prepend 'return_' to get the property name)",
                           dim=True, italic=True, colour='forest_green',
                           enabled=format)
            final_message += f"{title}\n{mess}\n"

        # Individual return flags
        mess = ''
        flags = {f'return_{ff}': vv for ff, vv in self._return_flags.items()
                 if f'return_{ff}' not in veto}
        for i, (flag, vv) in enumerate(flags.items()):
            prefix = "└" if i == len(flags) - 1 else "├"
            val = getattr(self, flag)
            name = f"{flag.replace('return_', '').replace('_', ' ')}:"
            mess += f"  {prefix} {name:15} {val}\n"
            subflags = [f'return_{ff}' for ff in vv
                        if f'return_{ff}' not in veto]
            for j, subflag in enumerate(subflags):
                subprefix = "└" if j == len(subflags) - 1 else "├"
                subval = getattr(self, subflag)
                subname = subflag.replace('return_', '')
                subname = subname.replace('other_', '').replace('_', ' ')
                subname = f"{subname}:"
                preprefix = " " if i == len(self._return_flags) - 1 else "│"
                mess += f"  {preprefix}   {subprefix} {subname:11} {subval}\n"
        if mess != '':
            title = "Individual return flags:"
            title = style(f"{title:25}", bold=True, colour='forest_green',
                          enabled=format)
            title += style(" (prepend 'return_' to get the property name)",
                           dim=True, italic=True, colour='forest_green',
                           enabled=format)
            final_message += f"{title}\n{mess}\n"

        # Physics processes
        mess = ''
        flags = [f'include_{pp}' for pp in self._include_processes
                 if f'include_{pp}' not in veto]
        for i, flag in enumerate(flags):
            prefix = "└" if i == len(flags) - 1 else "├"
            val = getattr(self, flag)
            name = f"{flag.replace('_', ' ')}:"
            mess += f"  {prefix} {name:32} {val}\n"
        if mess != '':
            title = "Physics processes:"
            title = style(title, bold=True, colour='forest_green',
                          enabled=format)
            title += style(" (prepend 'include_' to get the property name)",
                           dim=True, italic=True, colour='forest_green',
                           enabled=format)
            final_message += f"{title}\n{mess}\n"

        # Energy cuts
        mess = ''
        flags = [f'{ff}_cut' for ff in self._cut_definitions
                 if f'{ff}_cut' not in veto]
        for i, flag in enumerate(flags):
            prefix = "└" if i == len(flags) - 1 else "├"
            val = getattr(self, flag)
            name = f"{flag.replace('_cut', '').replace('_', ' ')}:"
            mess += f"  {prefix} {name:24} {val}\n"
        if mess != '':
            title = "Energy/momentum cuts [eV]:"
            title = style(f"{title:25}", bold=True, colour='forest_green',
                          enabled=format)
            title += style(" (append '_cut' to get the property name)",
                           dim=True, italic=True, colour='forest_green',
                           enabled=format)
            final_message += f"{title}\n{mess}\n"

        return final_message

    @property
    def _leaf_flags(self):
        nested_flags  = [f"return_{mm}" for mm in self._return_modifiers]
        nested_flags += [
            f"return_{kkk}"
            for kk, vv in self._return_flags.items()
            for kkk in (vv or [kk])
        ]
        nested_flags += [f"{cc}_cut" for cc in self._cut_definitions]
        nested_flags += [f"include_{pp}" for pp in self._include_processes]
        return nested_flags

    def _get_raw_settings(self):
        return {
            name: getattr(self, f"_{name}")
            for name in self._leaf_flags
        }

    def _set_raw_settings(self, settings):
        for name, value in settings.items():
            setattr(self, name, value)


    # ==========================
    # === Reference Particle ===
    # =========================

    @property
    def particle_ref(self):
        if self._engine.particle_ref is None:
            return xt.Particles()
        return self._engine.particle_ref

    @property
    def ref_id(self):
        return self.particle_ref.pdg_id[0]

    @property
    def ref_p0c(self):
        return self.particle_ref.p0c[0]

    @property
    def ref_is_lepton(self):
        return pdg.is_lepton(self.ref_id)

    @property
    def ref_is_proton(self):
        return pdg.is_proton(self.ref_id)

    @property
    def ref_is_ion(self):
        return pdg.is_ion(self.ref_id)


    # =====================
    # === Return groups ===
    # =====================

    @property
    def return_all(self):
        return (
            self.return_neutral and
            self.return_photons and
            self.return_electrons and
            self.return_muons and
            self.return_tauons and
            self.return_neutrinos and
            self.return_protons and
            self.return_neutrons and
            self.return_other_baryons and
            self.return_pions and
            self.return_kaons and
            self.return_other_mesons and
            self.return_ions
        )

    @return_all.setter
    def return_all(self, val):
        if val is True:
            self.return_neutral = True
            self.return_photons = True
            self.return_electrons = True
            self.return_muons = True
            self.return_tauons = True
            self.return_neutrinos = True
            self.return_protons = True
            self.return_neutrons = True
            self.return_other_baryons = True
            self.return_pions = True
            self.return_kaons = True
            self.return_other_mesons = True
            self.return_ions = True
        elif val is None or val is False:
            # Default settings
            self.return_neutral = None
            self.return_photons = None
            self.return_electrons = None
            self.return_muons = None
            self.return_tauons = None
            self.return_neutrinos = None
            self.return_protons = None
            self.return_neutrons = None
            self.return_other_baryons = None
            self.return_pions = None
            self.return_kaons = None
            self.return_other_mesons = None
            self.return_ions = None
        else:
            self._error("`return_all` has to be a boolean!")

    @property
    def return_all_charged(self):
        return (
            not self.return_neutral and
            self.return_electrons and
            self.return_muons and
            self.return_tauons and
            self.return_protons and
            self.return_other_baryons and
            self.return_pions and
            self.return_kaons and
            self.return_other_mesons and
            self.return_ions
        )

    @return_all_charged.setter
    def return_all_charged(self, val):
        if val is True:
            self.return_neutral = False
            self.return_photons = False
            self.return_electrons = True
            self.return_muons = True
            self.return_tauons = True
            self.return_neutrinos = False
            self.return_protons = True
            self.return_neutrons = False
            self.return_other_baryons = True
            self.return_pions = True
            self.return_kaons = True
            self.return_other_mesons = True
            self.return_ions = True
        elif val is not None and val is not False:
            self._error("`return_all_charged` has to be a boolean!")

    @property
    def return_none(self):
        return not (
            self.return_neutral or
            self.return_photons or
            self.return_electrons or
            self.return_muons or
            self.return_tauons or
            self.return_neutrinos or
            self.return_protons or
            self.return_neutrons or
            self.return_other_baryons or
            self.return_pions or
            self.return_kaons or
            self.return_other_mesons or
            self.return_ions
        )

    @return_none.setter
    def return_none(self, val):
        if val is True:
            self.return_neutral = False
            self.return_photons = False
            self.return_electrons = False
            self.return_muons = False
            self.return_tauons = False
            self.return_neutrinos = False
            self.return_protons = False
            self.return_neutrons = False
            self.return_other_baryons = False
            self.return_pions = False
            self.return_kaons = False
            self.return_other_mesons = False
            self.return_ions = False
        elif val is None or val is False:
            # Default settings
            self.return_neutral = None
            self.return_photons = None
            self.return_electrons = None
            self.return_muons = None
            self.return_tauons = None
            self.return_neutrinos = None
            self.return_protons = None
            self.return_neutrons = None
            self.return_other_baryons = None
            self.return_pions = None
            self.return_kaons = None
            self.return_other_mesons = None
            self.return_ions = None
        else:
            self._error("`return_none` has to be a boolean!")

    @property
    def return_neutral(self):
        if self._return_neutral is None:
            return False
        return self._return_neutral

    @return_neutral.setter
    def return_neutral(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_neutral` has to be a boolean!")
        self._return_neutral = val

    @property
    def return_leptons(self):
        ret = [getattr(self, f'return_{attr}')
               for attr in self._return_flags['leptons']]
        if all(ret):
            return True
        elif not any(ret):
            return False
        else:
            return None

    @return_leptons.setter
    def return_leptons(self, val):
        if not isinstance(val, bool) and val is not None:
            self._error("`return_leptons` has to be a boolean!")
        for attr in self._return_flags['leptons']:
            setattr(self, f'return_{attr}', val)

    @property
    def return_mesons(self):
        ret = [getattr(self, f'return_{attr}')
               for attr in self._return_flags['mesons']]
        if all(ret):
            return True
        elif not any(ret):
            return False
        else:
            return None

    @return_mesons.setter
    def return_mesons(self, val):
        if not isinstance(val, bool) and val is not None:
            self._error("`return_mesons` has to be a boolean!")
        for attr in self._return_flags['mesons']:
            setattr(self, f'return_{attr}', val)

    @property
    def return_baryons(self):
        ret = [getattr(self, f'return_{attr}')
               for attr in self._return_flags['baryons']]
        if all(ret):
            return True
        elif not any(ret):
            return False
        else:
            return None

    @return_baryons.setter
    def return_baryons(self, val):
        if not isinstance(val, bool) and val is not None:
            self._error("`return_baryons` has to be a boolean!")
        for attr in self._return_flags['baryons']:
            setattr(self, f'return_{attr}', val)


    # ==========================
    # === Return individuals ===
    # ==========================

    @property
    def return_photons(self):
        if self._return_photons is None:
            return self.return_neutral
        return self._return_photons

    @return_photons.setter
    def return_photons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_photons` has to be a boolean!")
        self._return_photons = val

    @property
    def return_electrons(self):
        if self._return_electrons is None:
            return self.ref_is_lepton
        return self._return_electrons

    @return_electrons.setter
    def return_electrons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_electrons` has to be a boolean!")
        self._return_electrons = val

    @property
    def return_muons(self):
        if self._return_muons is None:
            return self.ref_is_lepton
        return self._return_muons

    @return_muons.setter
    def return_muons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_muons` has to be a boolean!")
        self._return_muons = val

    @property
    def return_tauons(self):
        if self._return_tauons is None:
            return self.ref_is_lepton
        return self._return_tauons

    @return_tauons.setter
    def return_tauons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_tauons` has to be a boolean!")
        self._return_tauons = val

    @property
    def return_neutrinos(self):
        if self._return_neutrinos is None:
            return self.ref_is_lepton and self.return_neutral
        return self._return_neutrinos

    @return_neutrinos.setter
    def return_neutrinos(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_neutrinos` has to be a boolean!")
        self._return_neutrinos = val

    @property
    def return_protons(self):
        if self._return_protons is None:
            return self.ref_is_proton or self.ref_is_ion
        return self._return_protons

    @return_protons.setter
    def return_protons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_protons` has to be a boolean!")
        self._return_protons = val

    @property
    def return_neutrons(self):
        if self._return_neutrons is None:
            res = self.ref_is_proton or self.ref_is_ion
            return res and self.return_neutral
        return self._return_neutrons

    @return_neutrons.setter
    def return_neutrons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_neutrons` has to be a boolean!")
        self._return_neutrons = val

    @property
    def return_other_baryons(self):
        if self._return_other_baryons is None:
            return False
        return self._return_other_baryons

    @return_other_baryons.setter
    def return_other_baryons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_other_baryons` has to be a boolean!")
        self._return_other_baryons = val

    @property
    def return_pions(self):
        if self._return_pions is None:
            return False
        return self._return_pions

    @return_pions.setter
    def return_pions(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_pions` has to be a boolean!")
        self._return_pions = val

    @property
    def return_kaons(self):
        if self._return_kaons is None:
            return False
        return self._return_kaons

    @return_kaons.setter
    def return_kaons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_kaons` has to be a boolean!")
        self._return_kaons = val

    @property
    def return_other_mesons(self):
        if self._return_other_mesons is None:
            return False
        return self._return_other_mesons

    @return_other_mesons.setter
    def return_other_mesons(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_other_mesons` has to be a boolean!")
        self._return_other_mesons = val

    @property
    def return_ions(self):
        if self._return_ions is None:
            return self.ref_is_ion
        return self._return_ions

    @return_ions.setter
    def return_ions(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`return_ions` has to be a boolean!")
        self._return_ions = val


    # =====================
    # === Momentum cuts ===
    # =====================

    @property
    def hadron_lower_momentum_cut(self):
        if self._hadron_lower_momentum_cut is None:
            val = self.ref_p0c / 10
            if self.ref_is_ion:
                _, A, _, _ = pdg.get_properties_from_pdg_id(self.ref_id)
                val /= A
            return val
        return self._hadron_lower_momentum_cut

    @hadron_lower_momentum_cut.setter
    def hadron_lower_momentum_cut(self, val):
        if val is not None:
            if not isinstance(val, Number) or val < 0:
                self._error("`hadron_lower_momentum_cut` has to be a "
                            "non-negative number!")
            elif val < 1.e9:
                self._print(
                    f"Warning: Hadron lower momentum cut of {val/1.e9}GeV "
                    "is very low and will result in very long computation "
                    "times."
                )
        self._hadron_lower_momentum_cut = val

    @property
    def photon_lower_momentum_cut(self):
        if self._photon_lower_momentum_cut is None:
            return self.ref_p0c / 1000
        return self._photon_lower_momentum_cut

    @photon_lower_momentum_cut.setter
    def photon_lower_momentum_cut(self, val):
        if val is not None:
            if not isinstance(val, Number) or val < 0:
                self._error("`photon_lower_momentum_cut` has to be a "
                            "non-negative number!")
            elif val < 1.e3:
                self._print(
                    f"Warning: Photon lower momentum cut of {val/1.e3}keV "
                    "is very low and will result in very long computation "
                    "times."
                )
        self._photon_lower_momentum_cut = val

    @property
    def electron_lower_momentum_cut(self):
        if self._electron_lower_momentum_cut is None:
            if self.ref_is_lepton:
                return self.ref_p0c / 10
            else:
                return self.ref_p0c / 1000
        return self._electron_lower_momentum_cut

    @electron_lower_momentum_cut.setter
    def electron_lower_momentum_cut(self, val):
        if val is not None:
            if not isinstance(val, Number) or val < 0:
                self._error("`electron_lower_momentum_cut` has to be a "
                            "non-negative number!")
            elif val < 1.e6:
                self._print(
                    f"Warning: Electron lower momentum cut of {val/1.e6}MeV "
                    "is very low and will result in very long computation "
                    "times."
                )
        self._electron_lower_momentum_cut = val

    @property
    def relative_energy_cut(self):
        if self._relative_energy_cut is None:
            return 0.1
        return self._relative_energy_cut

    @relative_energy_cut.setter
    def relative_energy_cut(self, val):
        if val is not None:
            if not isinstance(val, Number) or val < 0:
                self._error("`relative_energy_cut` has to be a "
                            "non-negative number!")
            elif val < 1.e6:
                self._print(
                    f"Warning: Relative energy cut of {val/1.e6}MeV "
                    "is very low and will result in very long computation "
                    "times."
                )
        self._relative_energy_cut = val


    # ========================
    # === Physics settings ===
    # ========================

    @property
    def include_showers(self):
        if self._include_showers is None:
            return self.ref_is_lepton or self.return_all
        return self._include_showers

    @include_showers.setter
    def include_showers(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`include_showers` has to be a boolean!")
        self._include_showers = val

    @property
    def include_single_coulomb(self):
        if self._include_single_coulomb is None:
            return True
        return self._include_single_coulomb

    @include_single_coulomb.setter
    def include_single_coulomb(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`include_single_coulomb` has to be a boolean!")
        self._include_single_coulomb = val

    @property
    def include_multiple_coulomb(self):
        if self._include_multiple_coulomb is None:
            return True
        return self._include_multiple_coulomb

    @include_multiple_coulomb.setter
    def include_multiple_coulomb(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`include_multiple_coulomb` has to be a boolean!")
        self._include_multiple_coulomb = val

    @property
    def include_ionisation_fluctuations(self):
        if self._include_ionisation_fluctuations is None:
            return True
        return self._include_ionisation_fluctuations

    @include_ionisation_fluctuations.setter
    def include_ionisation_fluctuations(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`include_ionisation_fluctuations` has to be a boolean!")
        self._include_ionisation_fluctuations = val

    @property
    def include_pair_production(self):
        if self._include_pair_production is None:
            return True
        return self._include_pair_production

    @include_pair_production.setter
    def include_pair_production(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`include_pair_production` has to be a boolean!")
        self._include_pair_production = val

    @property
    def include_bremsstrahlung(self):
        if self._include_bremsstrahlung is None:
            return True
        return self._include_bremsstrahlung

    @include_bremsstrahlung.setter
    def include_bremsstrahlung(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`include_bremsstrahlung` has to be a boolean!")
        self._include_bremsstrahlung = val

    @property
    def include_elastic(self):
        if self._include_elastic is None:
            return True
        return self._include_elastic

    @include_elastic.setter
    def include_elastic(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`include_elastic` has to be a boolean!")
        self._include_elastic = val

    @property
    def include_inelastic(self):
        if self._include_inelastic is None:
            return True
        return self._include_inelastic

    @include_inelastic.setter
    def include_inelastic(self, val):
        if val is not None and not isinstance(val, bool):
            self._error("`include_inelastic` has to be a boolean!")
        self._include_inelastic = val


    # ======================
    # === Public Methods ===
    # ======================

    def reset(self):
        # Set flags to default
        for flag in [
            'return_all',
            *[f'{ff}_cut' for ff in self._cut_definitions],
            *[f'include_{pp}' for pp in self._include_processes],
        ]:
            if flag not in self._engine._physics_settings_veto_list:
                setattr(self, flag, None)


    def mask_particle_return_types(self, pdg_id, q_new):
        if self.return_all:
            # Allow everything and exclude
            mask_new = np.ones_like(pdg_id, dtype=bool)
        else:
            # Allow nothing and include
            mask_new = np.zeros_like(pdg_id, dtype=bool)

        # General categories
        is_ion = np.abs(pdg_id) > 1000000000
        mask_new[is_ion] = self.return_ions

        # PDG ID of mesons: from .*0XX. where X != 0 and . is any digit
        # restrict to |pdg_id| < 1e9 to not clash with ions
        mask_new[~is_ion & (pdg_id > 0) & (pdg_id // 10 % 10 != 0)
                 & (pdg_id // 100 % 10 != 0) & (pdg_id // 1000 % 10 == 0)
                ] = self.return_other_mesons
        mask_new[~is_ion & (pdg_id < 0) & (-pdg_id // 10 % 10 != 0)
                 & (-pdg_id // 100 % 10 != 0) & (-pdg_id // 1000 % 10 == 0)
                ] = self.return_other_mesons

        # PDG ID of baryons: from XXX. where X != 0 and . is any digit
        # restrict to |pdg_id| < 1e9 to not clash with ions
        mask_new[~is_ion & (pdg_id > 1000) & (pdg_id < 9000)
                 & (pdg_id // 10 % 10 != 0) & (pdg_id // 100 % 10 != 0)
                 & (pdg_id // 1000 % 10 != 0)
                ] = self.return_other_baryons    # PDG ID of from XX0X is a diquark
        mask_new[~is_ion & (pdg_id < -1000) & (pdg_id > -9000)
                 & (-pdg_id // 10 % 10 != 0) & (-pdg_id // 100 % 10 != 0)
                 & (-pdg_id // 1000 % 10 != 0)
                ] = self.return_other_baryons

        if not self.return_neutral:
            # General modifier, has to be before more specific return types,
            # as other neutral particles might have been specifically activated.
            mask_new[np.abs(q_new) < 1.e-12] = False

        mask_new[np.abs(pdg_id) == 22] = self.return_photons
        mask_new[np.abs(pdg_id) == 11] = self.return_electrons
        mask_new[np.abs(pdg_id) == 12] = self.return_electrons and self.return_neutrinos
        mask_new[np.abs(pdg_id) == 13] = self.return_muons
        mask_new[np.abs(pdg_id) == 14] = self.return_muons and self.return_neutrinos
        mask_new[np.abs(pdg_id) == 15] = self.return_tauons
        mask_new[np.abs(pdg_id) == 16] = self.return_tauons and self.return_neutrinos
        mask_new[np.abs(pdg_id) == 211] = self.return_pions
        mask_new[np.abs(pdg_id) == 111] = self.return_pions and self.return_neutral
        mask_new[np.abs(pdg_id) == 321] = self.return_kaons
        mask_new[np.abs(pdg_id) == 130] = self.return_kaons and self.return_neutral
        mask_new[np.abs(pdg_id) == 310] = self.return_kaons and self.return_neutral
        mask_new[np.abs(pdg_id) == 311] = self.return_kaons and self.return_neutral
        mask_new[np.abs(pdg_id) == 2212] = self.return_protons
        mask_new[np.abs(pdg_id) == 2112] = self.return_neutrons

        return mask_new


    # ======================
    # === Private Methods ===
    # ======================

    def __getattribute__(self, item):
        # Always use base lookup inside this method
        obj_get = object.__getattribute__
        try:
            engine = obj_get(self, "_engine")
        except AttributeError:
            # _engine does not exist yet during construction
            engine = None

        if item.startswith("return_"):
            # compute flags WITHOUT going through self.<property>
            global_flags     = obj_get(self, "_global_return_flags")
            return_flags     = obj_get(self, "_return_flags")
            return_modifiers = obj_get(self, "_return_modifiers")

            all_flags = list(global_flags)
            all_flags += list(return_flags.keys())
            all_flags += list(return_modifiers)
            all_flags += [fff for ff in return_flags.values() for fff in ff]

            if item not in {f"return_{ff}" for ff in all_flags}:
                if engine is not None:
                    engine.stop()
                raise AttributeError(f"Return flag '{item}' does not exist!")

        elif not item.startswith('_') and item.endswith("_cut"):
            cut_defs = obj_get(self, "_cut_definitions")
            if item not in {f"{ff}_cut" for ff in cut_defs}:
                if engine is not None:
                    engine.stop()
                raise AttributeError(f"Cut definition '{item}' does not exist!")

        elif item.startswith("include_"):
            include_processes = obj_get(self, "_include_processes")
            if item not in {f"include_{pp}" for pp in include_processes}:
                if engine is not None:
                    engine.stop()
                raise AttributeError(f"Physics process '{item}' does not exist!")

        if engine is not None and item in engine._physics_settings_veto_list:
            engine.stop()
            raise AttributeError(
                f"{engine.name.capitalize()} does not support "
                f"physics setting '{item}'"
            )

        return obj_get(self, item)


    def __setattr__(self, name, value):
        engine = self._engine

        if name in engine._physics_settings_veto_list:
            engine.stop()
            raise AttributeError(
                f"{engine.name.capitalize()} does not support "
                f"physics setting '{name}'"
            )

        if not hasattr(self, name):
            engine.stop()
            raise AttributeError(
                f"{self.__class__.__name__} object has no "
                f"attribute '{name}'"
            )

        super().__setattr__(name, value)

    @classmethod
    @contextmanager
    def _in_constructor(cls, self=None):
        original_setattr = cls.__setattr__
        def new_setattr(self, *args, **kwargs):
            return super().__setattr__( *args, **kwargs)
        cls.__setattr__ = new_setattr
        try:
            yield
        finally:
            cls.__setattr__ = original_setattr
