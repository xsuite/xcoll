# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

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

    def __init__(self, engine):
        with self.__class__._in_constructor(self):
            self._engine = engine
            self.reset()

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))} (use .show() " \
             + f"to see the contents)>"

    def __str__(self):
        return self._str(format=False)

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

    def update(self):
        # Update settings when they are set to use defaults
        for kk, vv in self.__dict__.items():
            if kk.startswith('_') and not kk.startswith('__') and kk.endswith('_use_default'):
                if vv:
                    prop_name = kk[1:-12]
                    if prop_name not in self._engine._physics_settings_veto_list:
                        setattr(self, prop_name, None)

    def reset(self):
        # Set flags to default
        for flag in [
            'return_all',
            *[f'{ff}_cut' for ff in self._cut_definitions],
            *[f'include_{pp}' for pp in self._include_processes],
        ]:
            if flag not in self._engine._physics_settings_veto_list:
                setattr(self, flag, None)

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
                subname = f"{subflag.replace('return_', '').replace('other_', '').replace('_', ' ')}:"
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
        return (self._return_neutral and
                self._return_photons and
                self._return_electrons and
                self._return_muons and
                self._return_tauons and
                self._return_neutrinos and
                self._return_protons and
                self._return_neutrons and
                self._return_other_baryons and
                self._return_pions and
                self._return_kaons and
                self._return_other_mesons and
                self._return_ions)

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
            self._engine.stop()
            raise ValueError("`return_all` has to be a boolean!")

    @property
    def return_all_charged(self):
        return (not self._return_neutral and
                self._return_electrons and
                self._return_muons and
                self._return_tauons and
                self._return_protons and
                self._return_other_baryons and
                self._return_pions and
                self._return_kaons and
                self._return_other_mesons and
                self._return_ions)

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
            self._engine.stop()
            raise ValueError("`return_all_charged` has to be a boolean!")

    @property
    def return_none(self):
        return not (self._return_neutral or
                self._return_photons or
                self._return_electrons or
                self._return_muons or
                self._return_tauons or
                self._return_neutrinos or
                self._return_protons or
                self._return_neutrons or
                self._return_other_baryons or
                self._return_pions or
                self._return_kaons or
                self._return_other_mesons or
                self._return_ions)

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
            self._engine.stop()
            raise ValueError("`return_none` has to be a boolean!")

    @property
    def return_neutral(self):
        return self._return_neutral

    @return_neutral.setter
    def return_neutral(self, val):
        if val is None:
            val = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_neutral` has to be a boolean!")
        self._return_neutral = val
        # This is a modifier flag; need to update default dependent flags
        self.update()
        if val is False:
            self.return_photons = False
            self.return_neutrinos = False
            self.return_neutrons = False

    @property
    def return_leptons(self):
        ret = [getattr(self, f'return_{attr}') for attr in self._return_flags['leptons']]
        if all(ret):
            return True
        elif not any(ret):
            return False
        else:
            return None

    @return_leptons.setter
    def return_leptons(self, val):
        if not isinstance(val, bool) and val is not None:
            self._engine.stop()
            raise ValueError("`return_leptons` has to be a boolean!")
        for attr in self._return_flags['leptons']:
            setattr(self, f'return_{attr}', val)

    @property
    def return_mesons(self):
        ret = [getattr(self, f'return_{attr}') for attr in self._return_flags['mesons']]
        if all(ret):
            return True
        elif not any(ret):
            return False
        else:
            return None

    @return_mesons.setter
    def return_mesons(self, val):
        if not isinstance(val, bool) and val is not None:
            self._engine.stop()
            raise ValueError("`return_mesons` has to be a boolean!")
        for attr in self._return_flags['mesons']:
            setattr(self, f'return_{attr}', val)

    @property
    def return_baryons(self):
        ret = [getattr(self, f'return_{attr}') for attr in self._return_flags['baryons']]
        if all(ret):
            return True
        elif not any(ret):
            return False
        else:
            return None

    @return_baryons.setter
    def return_baryons(self, val):
        if not isinstance(val, bool) and val is not None:
            self._engine.stop()
            raise ValueError("`return_baryons` has to be a boolean!")
        for attr in self._return_flags['baryons']:
            setattr(self, f'return_{attr}', val)


    # ==========================
    # === Return individuals ===
    # ==========================

    @property
    def return_photons(self):
        return self._return_photons

    @return_photons.setter
    def return_photons(self, val):
        if val is None:
            self._return_photons_use_default = True
            val = self.return_neutral
        else:
            self._return_photons_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_photons` has to be a boolean!")
        self._return_photons = val

    @property
    def return_electrons(self):
        return self._return_electrons

    @return_electrons.setter
    def return_electrons(self, val):
        if val is None:
            self._return_electrons_use_default = True
            val = self.ref_is_lepton
        else:
            self._return_electrons_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_electrons` has to be a boolean!")
        self._return_electrons = val

    @property
    def return_muons(self):
        return self._return_muons

    @return_muons.setter
    def return_muons(self, val):
        if val is None:
            self._return_muons_use_default = True
            val = self.ref_is_lepton
        else:
            self._return_muons_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_muons` has to be a boolean!")
        self._return_muons = val

    @property
    def return_tauons(self):
        return self._return_tauons

    @return_tauons.setter
    def return_tauons(self, val):
        if val is None:
            self._return_tauons_use_default = True
            val = self.ref_is_lepton
        else:
            self._return_tauons_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_tauons` has to be a boolean!")
        self._return_tauons = val

    @property
    def return_neutrinos(self):
        return self._return_neutrinos

    @return_neutrinos.setter
    def return_neutrinos(self, val):
        if val is None:
            self._return_neutrinos_use_default = True
            val = self.ref_is_lepton and self.return_neutral
        else:
            self._return_neutrinos_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_neutrinos` has to be a boolean!")
        self._return_neutrinos = val

    @property
    def return_protons(self):
        return self._return_protons

    @return_protons.setter
    def return_protons(self, val):
        if val is None:
            self._return_protons_use_default = True
            val = self.ref_is_proton or self.ref_is_ion
        else:
            self._return_protons_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_protons` has to be a boolean!")
        self._return_protons = val

    @property
    def return_neutrons(self):
        return self._return_neutrons

    @return_neutrons.setter
    def return_neutrons(self, val):
        if val is None:
            self._return_neutrons_use_default = True
            val = (self.ref_is_proton or self.ref_is_ion) and self.return_neutral
        else:
            self._return_neutrons_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_neutrons` has to be a boolean!")
        self._return_neutrons = val

    @property
    def return_other_baryons(self):
        return self._return_other_baryons

    @return_other_baryons.setter
    def return_other_baryons(self, val):
        if val is None:
            val = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_other_baryons` has to be a boolean!")
        self._return_other_baryons = val

    @property
    def return_pions(self):
        return self._return_pions

    @return_pions.setter
    def return_pions(self, val):
        if val is None:
            val = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_pions` has to be a boolean!")
        self._return_pions = val

    @property
    def return_kaons(self):
        return self._return_kaons

    @return_kaons.setter
    def return_kaons(self, val):
        if val is None:
            val = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_kaons` has to be a boolean!")
        self._return_kaons = val

    @property
    def return_other_mesons(self):
        return self._return_other_mesons

    @return_other_mesons.setter
    def return_other_mesons(self, val):
        if val is None:
            val = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_other_mesons` has to be a boolean!")
        self._return_other_mesons = val

    @property
    def return_ions(self):
        return self._return_ions

    @return_ions.setter
    def return_ions(self, val):
        if val is None:
            self._return_ions_use_default = True
            val = self.ref_is_ion
        else:
            self._return_ions_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`return_ions` has to be a boolean!")
        self._return_ions = val


    # =====================
    # === Momentum cuts ===
    # =====================

    @property
    def hadron_lower_momentum_cut(self):
        return self._hadron_lower_momentum_cut

    @hadron_lower_momentum_cut.setter
    def hadron_lower_momentum_cut(self, val):
        self._hadron_lower_momentum_cut_use_default = False
        if val is None:
            self._hadron_lower_momentum_cut_use_default = True
            val = self.ref_p0c / 10
            if self.ref_is_ion:
                _, A, _, _ = pdg.get_properties_from_pdg_id(self.ref_id)
                val /= A
        elif not isinstance(val, Number) or val < 0:
            self._engine.stop()
            raise ValueError("`hadron_lower_momentum_cut` has to be a non-negative number!")
        elif val < 1.e9:
            self._print(f"Warning: Hadron lower momentum cut of {val/1.e9}GeV "
                       + "is very low and will result in very long computation times.")
        self._hadron_lower_momentum_cut = val

    @property
    def photon_lower_momentum_cut(self):
        return self._photon_lower_momentum_cut

    @photon_lower_momentum_cut.setter
    def photon_lower_momentum_cut(self, val):
        self._photon_lower_momentum_cut_use_default = False
        if val is None:
            self._photon_lower_momentum_cut_use_default = True
            val = self.ref_p0c * 1e-3
        elif not isinstance(val, Number) or val < 0:
            self._engine.stop()
            raise ValueError("`photon_lower_momentum_cut` has to be a non-negative number!")
        elif val < 1.e3:
            self._print(f"Warning: Photon lower momentum cut of {val/1.e3}keV "
                       + "is very low and will result in very long computation times.")
        self._photon_lower_momentum_cut = val

    @property
    def electron_lower_momentum_cut(self):
        return self._electron_lower_momentum_cut

    @electron_lower_momentum_cut.setter
    def electron_lower_momentum_cut(self, val):
        self._electron_lower_momentum_cut_use_default = False
        if val is None:
            self._electron_lower_momentum_cut_use_default = True
            if self.ref_is_lepton:
                val = self.ref_p0c / 10
            else:
                val = self.ref_p0c * 1e-3
        elif not isinstance(val, Number) or val < 0:
            self._engine.stop()
            raise ValueError("`electron_lower_momentum_cut` has to be a non-negative number!")
        elif val < 1.e6:
            self._print(f"Warning: Electron lower momentum cut of {val/1.e6}MeV "
                       + "is very low and will result in very long computation times.")
        self._electron_lower_momentum_cut = val

    @property
    def relative_energy_cut(self):
        return self._relative_energy_cut

    @relative_energy_cut.setter
    def relative_energy_cut(self, val):
        self._relative_energy_cut_use_default = False
        if val is None:
            self._relative_energy_cut_use_default = True
            val = 0.1
        elif not isinstance(val, Number) or val <= 0:
            self._engine.stop()
            raise ValueError("`relative_energy_cut` has to be a strictly positive number!")
        elif val < 1e-6:
            self._print(f"Warning: Relative energy cut of {val} is very low and will "
                       + "result in very long computation times.")
        self._relative_energy_cut = val


    # ========================
    # === Physics settings ===
    # ========================

    @property
    def include_showers(self):
        return self._include_showers

    @include_showers.setter
    def include_showers(self, val):
        if val is None:
            self._include_showers_use_default = True
            val = True if self.ref_is_lepton else False
            val = True if self.return_all else val
        else:
            self._include_showers_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`include_showers` has to be a boolean!")
        self._include_showers = val

    @property
    def include_single_coulomb(self):
        return self._include_single_coulomb

    @include_single_coulomb.setter
    def include_single_coulomb(self, val):
        if val is None:
            self._include_single_coulomb_use_default = True
            val = True
        else:
            self._include_single_coulomb_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`include_single_coulomb` has to be a boolean!")
        self._include_single_coulomb = val

    @property
    def include_multiple_coulomb(self):
        return self._include_multiple_coulomb

    @include_multiple_coulomb.setter
    def include_multiple_coulomb(self, val):
        if val is None:
            self._include_multiple_coulomb_use_default = True
            val = True
        else:
            self._include_multiple_coulomb_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`include_multiple_coulomb` has to be a boolean!")
        self._include_multiple_coulomb = val

    @property
    def include_ionisation_fluctuations(self):
        return self._include_ionisation_fluctuations

    @include_ionisation_fluctuations.setter
    def include_ionisation_fluctuations(self, val):
        if val is None:
            self._include_ionisation_fluctuations_use_default = True
            val = True
        else:
            self._include_ionisation_fluctuations_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`include_ionisation_fluctuations` has to be a boolean!")
        self._include_ionisation_fluctuations = val

    @property
    def include_pair_production(self):
        return self._include_pair_production

    @include_pair_production.setter
    def include_pair_production(self, val):
        if val is None:
            self._include_pair_production_use_default = True
            val = True
        else:
            self._include_pair_production_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`include_pair_production` has to be a boolean!")
        self._include_pair_production = val

    @property
    def include_bremsstrahlung(self):
        return self._include_bremsstrahlung

    @include_bremsstrahlung.setter
    def include_bremsstrahlung(self, val):
        if val is None:
            self._include_bremsstrahlung_use_default = True
            val = True
        else:
            self._include_bremsstrahlung_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`include_bremsstrahlung` has to be a boolean!")
        self._include_bremsstrahlung = val

    @property
    def include_elastic(self):
        return self._include_elastic

    @include_elastic.setter
    def include_elastic(self, val):
        if val is None:
            self._include_elastic_use_default = True
            val = True
        else:
            self._include_elastic_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`include_elastic` has to be a boolean!")
        self._include_elastic = val

    @property
    def include_inelastic(self):
        return self._include_inelastic

    @include_inelastic.setter
    def include_inelastic(self, val):
        if val is None:
            self._include_inelastic_use_default = True
            val = True
        else:
            self._include_inelastic_use_default = False
        if not isinstance(val, bool):
            self._engine.stop()
            raise ValueError("`include_inelastic` has to be a boolean!")
        self._include_inelastic = val


    def __getattribute__(self, item):
        # Always use base lookup inside this method
        obj_get = object.__getattribute__

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
                engine = obj_get(self, "_engine")
                engine.stop()
                raise AttributeError(f"Return flag '{item}' does not exist!")

        elif not item.startswith('_') and item.endswith("_cut"):
            cut_defs = obj_get(self, "_cut_definitions")
            if item not in {f"{ff}_cut" for ff in cut_defs}:
                engine = obj_get(self, "_engine")
                engine.stop()
                raise AttributeError(f"Cut definition '{item}' does not exist!")

        elif item.startswith("include_"):
            include_processes = obj_get(self, "_include_processes")
            if item not in {f"include_{pp}" for pp in include_processes}:
                engine = obj_get(self, "_engine")
                engine.stop()
                raise AttributeError(f"Physics process '{item}' does not exist!")

        try:
            engine = obj_get(self, "_engine")
            if item in engine._physics_settings_veto_list:
                engine.stop()
                raise AttributeError(f"{engine.name.capitalize()} does not support "
                                     f"physics setting '{item}'")
        except AttributeError:
            pass

        return obj_get(self, item)

    def __setattr__(self, name, value):
        if not hasattr(self, name):
            self._engine.stop()
            raise AttributeError(f"{self.__class__.__name__} object has no "
                                 f"attribute '{name}'")
        if name in self._engine._physics_settings_veto_list:
            self._engine.stop()
            raise AttributeError(f"{self._engine.name.capitalize()} does not "
                                 f"support physics setting '{name}'")
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
