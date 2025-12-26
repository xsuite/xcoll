# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from numbers import Number

import xtrack as xt
import xtrack.particles.pdg as pdg


class PhysicsSettingsHelper:
    """Helper class to manage physics settings for scattering routines.
    """

    _return_flags = {'photons': [],
                     'leptons': ['electrons', 'muons', 'tauons', 'neutrinos'],
                     'baryons': ['protons', 'neutrons', 'other_baryons'],
                     'mesons':  ['pions', 'kaons', 'other_mesons'],
                     'ions':    []}
    _global_return_flags = ['all', 'all_charged', 'none']
    _return_modifiers = ['neutral']
    _cut_definitions = ['hadron_lower_momentum', 'photon_lower_momentum', 'electron_lower_momentum']

    def __init__(self, engine):
        self._engine = engine
        self._use_cuts = True
        # Set flags to default
        self.return_all = None
        self.hadron_lower_momentum_cut = None
        self.photon_lower_momentum_cut = None
        self.electron_lower_momentum_cut = None
        self.include_showers = None

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

    def update(self):
        # Update settings when they are set to use defaults
        for kk, vv in self.__dict__.items():
            if kk.startswith('_') and not kk.startswith('__') and kk.endswith('_use_default'):
                if vv:
                    prop_name = kk[1:-12]
                    setattr(self, prop_name, None)

    def show(self):
        print("Global return flags:")
        for i, ff in enumerate(self._global_return_flags):
            prefix = "└" if i == len(self._global_return_flags) - 1 else "├"
            val = getattr(self, f'return_{ff}')
            flag = f"{ff.replace('_', ' ')}:"
            print(f"  {prefix} {flag:15} {val}")
        print()
        print("Global return modifiers:")
        for i, ff in enumerate(self._return_modifiers):
            prefix = "└" if i == len(self._return_modifiers) - 1 else "├"
            val = getattr(self, f'return_{ff}')
            flag = f"{ff.replace('_', ' ')}:"
            print(f"  {prefix} {flag:15} {val}")
        print()
        print("Individual return flags:")
        for i, (ff, vv) in enumerate(self._return_flags.items()):
            prefix = "└" if i == len(self._return_flags) - 1 else "├"
            val = getattr(self, f'return_{ff}')
            flag = f"{ff.replace('_', ' ')}:"
            print(f"  {prefix} {flag:15} {val}")
            for j, subff in enumerate(vv):
                subprefix = "└" if j == len(vv) - 1 else "├"
                subval = getattr(self, f'return_{subff}')
                subflag = f"{subff.replace('other_', '').replace('_', ' ')}:"
                preprefix = " " if i == len(self._return_flags) - 1 else "│"
                print(f"  {preprefix}   {subprefix} {subflag:11} {subval}")
        if self._use_cuts:
            print()
            print("Energy cuts [eV]:")
            for i, ff in enumerate(self._cut_definitions):
                prefix = "└" if i == len(self._cut_definitions) - 1 else "├"
                val = getattr(self, f'{ff}_cut')
                flag = f"{ff.replace('_', ' ')}:"
                print(f"  {prefix} {flag:28} {val}")

    @property
    def all_return_flags(self):
        all_flags  = self._global_return_flags.copy()
        all_flags += self._return_flags.keys()
        all_flags += self._return_modifiers
        all_flags += [fff for ff in self._return_flags.values() for fff in ff]
        return [f'return_{ff}' for ff in all_flags]

    @property
    def all_cut_definitions(self):
        return [f'{ff}_cut' for ff in self._cut_definitions]

    def __getattribute__(self, item):
        # always use base lookup inside this method
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
                engine = obj_get(self, "_engine") # also fetch engine via base lookup to avoid re-entry surprises
                engine.stop()
                raise AttributeError(f"Return flag '{item}' does not exist!")

        elif not item.startswith('_') and item.endswith("_cut"):
            cut_defs = obj_get(self, "_cut_definitions")
            if item not in {f"{ff}_cut" for ff in cut_defs}:
                engine = obj_get(self, "_engine")
                engine.stop()
                raise AttributeError(f"Cut definition '{item}' does not exist!")

        return obj_get(self, item)


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
        if not self._use_cuts:
            return None
        return self._hadron_lower_momentum_cut

    @hadron_lower_momentum_cut.setter
    def hadron_lower_momentum_cut(self, val):
        if not self._use_cuts:
            self._hadron_lower_momentum_cut = None
            return
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
        if not self._use_cuts:
            return None
        return self._photon_lower_momentum_cut

    @photon_lower_momentum_cut.setter
    def photon_lower_momentum_cut(self, val):
        if not self._use_cuts:
            self._photon_lower_momentum_cut = None
            return
        self._photon_lower_momentum_cut_use_default = False
        if val is None:
            self._photon_lower_momentum_cut_use_default = True
            val = 1.e5
        elif not isinstance(val, Number) or val < 0:
            self._engine.stop()
            raise ValueError("`photon_lower_momentum_cut` has to be a non-negative number!")
        elif val < 1.e3:
            self._print(f"Warning: Photon lower momentum cut of {val/1.e3}keV "
                       + "is very low and will result in very long computation times.")
        self._photon_lower_momentum_cut = val

    @property
    def electron_lower_momentum_cut(self):
        if not self._use_cuts:
            return None
        return self._electron_lower_momentum_cut

    @electron_lower_momentum_cut.setter
    def electron_lower_momentum_cut(self, val):
        if not self._use_cuts:
            self._electron_lower_momentum_cut = None
            return
        self._electron_lower_momentum_cut_use_default = False
        if val is None:
            self._electron_lower_momentum_cut_use_default = True
            if self.ref_is_lepton:
                val = self.ref_p0c / 10
            else:
                val = 1.e6
        elif not isinstance(val, Number) or val < 0:
            self._engine.stop()
            raise ValueError("`electron_lower_momentum_cut` has to be a non-negative number!")
        elif val < 1.e6:
            self._print(f"Warning: Electron lower momentum cut of {val/1.e6}MeV "
                       + "is very low and will result in very long computation times.")
        self._electron_lower_momentum_cut = val

    @property
    def include_showers(self):
        if not self._use_cuts:
            return None
        return self._include_showers

    @include_showers.setter
    def include_showers(self, val):
        if not self._use_cuts:
            self._include_showers = None
            return
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
