# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from numbers import Number

import xtrack as xt
import xtrack.particles.pdg as pdg

from ..pretty_print import style
from ..xaux import track_construction, super_if_being_constructed


@track_construction
class PhysicsSettingsHelper:
    """Helper class to manage physics settings for scattering routines.
    """

    # List of all return flags, sorted by family.
    # Each family is a dictionary with the following possible keys:
    # - 'default': a function that takes the PhysicsSettingsHelper instance
    #               and returns a boolean indicating whether to return this
    #               family by default (depending on the reference particle).
    # - 'pdg_ids': a tuple of PDG IDs corresponding to this family.
    # - 'neutral_pdg_ids': a tuple of PDG IDs corresponding to neutral
    #                      particles in this family. Only use this if the
    #                      family contains both charged and neutral particles.
    # - 'selector': a function that takes a PDG ID and returns a boolean
    #               indicating whether this PDG ID belongs to this family.
    # The last three keys are used in the `mask_particle_return_types`` method
    # to determine whether a particle with a given PDG ID should be returned.
    #
    # The user can set/unset these flags by prepending it with 'return_'.

    _global_return_flags = ['all', 'all_charged', 'none']
    _return_flags = {
        'photons': {
            'photons': {
                'default': lambda self: self.ref_is_photon,
                'pdg_ids': (22,),
            },
        },
        'leptons': {
            'electrons': {
                'default': lambda self: self.ref_is_lepton,
                'pdg_ids': (-11, 11),
            },
            'muons': {
                'default': lambda self: self.ref_is_lepton,
                'pdg_ids': (-13, 13),
            },
            'tauons': {
                'default': lambda self: self.ref_is_lepton,
                'pdg_ids': (-15, 15),
            },
            'neutrinos': {
                'default': lambda self:
                                self.ref_is_lepton and self.return_neutral,
                'pdg_ids': (-16, -14, -12, 12, 14, 16),
            },
        },
        'baryons': {
            'protons': {
                'default': lambda self: self.ref_is_baryon or self.ref_is_ion,
                'pdg_ids': (-2212, 2212),
            },
            'neutrons': {
                'default': lambda self:
                                self.return_protons and self.return_neutral,
                'pdg_ids': (-2112, 2112),
            },
            'other_baryons': {
                'default': lambda self:
                                self.ref_is_baryon and not self.ref_is_nucleon,
                'selector': lambda pdg_id: (
                                pdg.is_baryon(pdg_id)
                                & ~np.isin(np.abs(pdg_id), [2212, 2112])
                            ),
            },
        },
        'mesons': {
            'pions': {
                'default': lambda self: self.ref_is_meson,
                'pdg_ids': (-211, 211),
                'neutral_pdg_ids': (111,),
            },
            'kaons': {
                'default': lambda self:
                                self.ref_is_meson and not self.ref_is_pion,
                'pdg_ids': (-321, 321),
                'neutral_pdg_ids': (-311, 130, 310, 311),
            },
            'other_mesons': {
                'default': lambda self:
                                self.ref_is_meson
                                and not self.ref_is_pion
                                and not self.ref_is_kaon,
                'selector': lambda pdg_id: (
                                pdg.is_meson(pdg_id)
                                & ~pdg.is_pion(pdg_id)
                                & ~pdg.is_kaon(pdg_id)
                            ),
            },
        },
        'ions': {
            'ions': {
                'default': lambda self: self.ref_is_ion,
                'selector': pdg.is_ion,
            },
        },
    }
    _return_modifiers = {
        'neutral': lambda self: self.ref_is_neutral,
    }
    _neutral_only_return_flags = {
        "photons",
        "neutrinos",
        "neutrons",
    }

    # Append '_cut' to get the property name for cut definitions
    _cut_definitions = [
        'hadron_lower_momentum',
        'photon_lower_momentum',
        'electron_lower_momentum',
        'relative_energy'
    ]

    # Prepend 'include_' to get the property name for the physics processes
    _include_processes = {
        'showers': lambda self: self.ref_is_lepton or self.return_all,
        'single_coulomb': True,
        'multiple_coulomb': True,
        'ionisation_fluctuations': True,
        'pair_production': True,
        'bremsstrahlung': True,
        'elastic': True,
        'inelastic': True
    }

    def __init__(self, engine):
        self._engine = engine
        self.reset()

    @property
    def all_flags(self):
        all_flags  = self._global_return_flags.copy()
        all_flags += list(self._return_flags)
        all_flags += list(self._return_modifiers)
        all_flags += [
            flag
            for group, flags in self._return_flags.items()
            for flag in flags.keys()
            if flag != group
        ]
        result  = [f"return_{ff}" for ff in all_flags]
        result += [f"{ff}_cut" for ff in self._cut_definitions]
        result += [f"include_{pp}" for pp in self._include_processes]
        return result


    # ==========================
    # === Reference Particle ===
    # ==========================

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
    def ref_is_nucleon(self):
        res  = (np.abs(self.ref_id) == 2212)  # proton
        res |= (np.abs(self.ref_id) == 2112)  # neutron
        return res
    @property
    def ref_is_baryon(self):
        return pdg.is_baryon(self.ref_id)
    @property
    def ref_is_ion(self):
        return pdg.is_ion(self.ref_id)
    @property
    def ref_is_pion(self):
        return pdg.is_pion(self.ref_id)
    @property
    def ref_is_kaon(self):
        return pdg.is_kaon(self.ref_id)
    @property
    def ref_is_meson(self):
        return pdg.is_meson(self.ref_id)
    @property
    def ref_is_lepton(self):
        return pdg.is_lepton(self.ref_id)
    @property
    def ref_is_photon(self):
        return pdg.is_photon(self.ref_id)
    @property
    def ref_is_neutral(self):
        return abs(self.particle_ref.q0) < 1e-12


    # ========================================
    # === Return flags (non-autogenerated) ===
    # ========================================

    @property
    def return_all(self):
        return (
            self.return_neutral
            and all(
                getattr(self, flag)
                for flag in self._return_leaf_flags
            )
        )

    @return_all.setter
    def return_all(self, val):
        if val is False:
            val = None
        elif val is not None and val is not True:
            self._error("`return_all` has to be a boolean!")
        self.return_neutral = val
        for flag in self._return_leaf_flags:
            setattr(self, flag, val)

    @property
    def return_all_charged(self):
        return (
            not self.return_neutral
            and all(
                getattr(self, flag)
                for flag in self._charged_return_leaf_flags
            )
            and not any(
                getattr(self, flag)
                for flag in self._neutral_only_return_leaf_flags
            )
        )

    @return_all_charged.setter
    def return_all_charged(self, val):
        if val is False:
            val = None
        elif val is not None and val is not True:
            self._error("`return_all_charged` has to be a boolean!")
        if val is None:
            self.return_neutral = None
            for flag in self._return_leaf_flags:
                setattr(self, flag, None)
            return
        self.return_neutral = False
        for flag in self._charged_return_leaf_flags:
            setattr(self, flag, True)
        for flag in self._neutral_only_return_leaf_flags:
            setattr(self, flag, False)

    @property
    def return_none(self):
        return not (
            self.return_neutral
            or any(
                getattr(self, flag)
                for flag in self._return_leaf_flags
            )
        )

    @return_none.setter
    def return_none(self, val):
        if val is False:
            val = None
        elif val is True:
            val = False
        elif val is not None:
            self._error("`return_none` has to be a boolean!")
        self.return_neutral = val
        for flag in self._return_leaf_flags:
            setattr(self, flag, val)

    def return_pdg_id(self, pdg_id):
        if not hasattr(pdg_id, "__iter__") or isinstance(pdg_id, str):
            pdg_id = [pdg_id]
        pdg_id = set([int(pid) for pid in pdg_id])
        self._extra_pdg_ids_to_return.update(pdg_id)
        self._extra_pdg_ids_to_kill -= pdg_id

    def dont_return_pdg_id(self, pdg_id):
        if not hasattr(pdg_id, "__iter__") or isinstance(pdg_id, str):
            pdg_id = [pdg_id]
        pdg_id = set([int(pid) for pid in pdg_id])
        self._extra_pdg_ids_to_kill.update(pdg_id)
        self._extra_pdg_ids_to_return -= pdg_id


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
            if not isinstance(val, Number) or val <= 0:
                self._error("`relative_energy_cut` has to be a "
                            "strictly positive number!")
            elif val < 1.e-6:
                self._print(
                    f"Warning: Relative energy cut of {val} is very low and "
                    "will result in very long computation times."
                )
        self._relative_energy_cut = val


    # ======================
    # === Public Methods ===
    # ======================

    def show(self, full=False):
        """Print the physics settings."""
        print(self._str(format=True))

    def reset(self):
        # Set flags to default
        for flag in [
            'return_all',
            *[f'{ff}_cut' for ff in self._cut_definitions],
            *[f'include_{pp}' for pp in self._include_processes],
        ]:
            if flag not in self._engine._physics_settings_veto_list:
                setattr(self, flag, None)
        self._extra_pdg_ids_to_return = set()
        self._extra_pdg_ids_to_kill = set()


    def pdg_id_is_returned(self, pdg_id, q=None):
        scalar = np.ndim(pdg_id) == 0
        pdg_id = np.atleast_1d(np.asarray(pdg_id, dtype=np.int64))
        result = np.zeros(pdg_id.shape, dtype=bool)

        explicit_return = np.isin(pdg_id, list(self._extra_pdg_ids_to_return))
        explicit_kill = np.isin(pdg_id, list(self._extra_pdg_ids_to_kill))
        # Exact overrides have highest priority and do not require
        # Xtrack to know anything else about the PDG ID.
        result[explicit_return] = True
        result[explicit_kill] = False

        unresolved = ~(explicit_return| explicit_kill)
        if np.any(unresolved):
            pids_un = pdg_id[unresolved]
            if q is None:
                try:
                    q_un = pdg.get_properties_from_pdg_id(pids_un)[0]
                except ValueError:
                    raise ValueError(f"Could not get charge for some PDG IDs. "
                                     f"Please provide `q` explicitly.")
            else:
                q = np.atleast_1d(np.asarray(q))
                if q.shape != pdg_id.shape:
                    raise ValueError("`q` must have the same shape "
                                     "as `pdg_id`.")
                q_un = q[unresolved]
            result[unresolved] = self.mask_particle_return_types(pids_un, q_un)

        if scalar:
            return bool(result[0])
        return result


    def mask_particle_return_types(self, pdg_id, q_new):
        pdg_id = np.asarray(pdg_id, dtype=np.int64)
        q_new = np.asarray(q_new)

        if self.return_all:
            # Allow everything and exclude
            mask_new = np.ones_like(pdg_id, dtype=bool)
        else:
            # Allow nothing and include
            mask_new = np.zeros_like(pdg_id, dtype=bool)

        # Families of particles
        for flags in self._return_flags.values():
            for flag, spec in flags.items():
                selector = spec.get("selector", None)
                if selector is not None:
                    value = getattr(self, f"return_{flag}")
                    mask_new[selector(pdg_id)] = value

        # Neutral particles
        if not self.return_neutral:
            # General modifier, has to be before more specific return types,
            # as other neutral particles might have been specifically activated.
            mask_new[np.abs(q_new) < 1.e-12] = False

        # Individual return types
        for flags in self._return_flags.values():
            for flag, spec in flags.items():
                value = getattr(self, f"return_{flag}")
                ids = spec.get("pdg_ids", ())
                if ids:
                    mask_new[np.isin(pdg_id, ids)] = value

                neutral_ids = spec.get("neutral_pdg_ids", ())
                if neutral_ids:
                    mask_new[np.isin(pdg_id, neutral_ids)] = (
                        value and self.return_neutral
                    )

        # Neutrinos are dealt with by lepton family
        mask_new[np.abs(pdg_id) == 12] = (
            self.return_neutrinos and self.return_electrons
        )
        mask_new[np.abs(pdg_id) == 14] = (
            self.return_neutrinos and self.return_muons
        )
        mask_new[np.abs(pdg_id) == 16] = (
            self.return_neutrinos and self.return_tauons
        )

        # Specific requested PDG IDs to return
        for pp in self._extra_pdg_ids_to_return:
            mask_new[pdg_id == pp] = True
        for pp in self._extra_pdg_ids_to_kill:
            mask_new[pdg_id == pp] = False

        return mask_new


    # =======================
    # === Private Methods ===
    # =======================

    def __repr__(self):
        return f"<{self.__class__.__name__} at {hex(id(self))}>"

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

        # Extra PDG IDs to return/kill
        if self._extra_pdg_ids_to_return or self._extra_pdg_ids_to_kill:
            title = "PDG ID overrides:"
            title = style(f"{title:25}", bold=True, colour='forest_green',
                          enabled=format)
            final_message += f"{title}\n"

            if self._extra_pdg_ids_to_return:
                values = ", ".join(
                    _format_pdg_id(pp)
                    for pp in sorted(self._extra_pdg_ids_to_return)
                )
                prefix = "├" if self._extra_pdg_ids_to_kill else "└"
                final_message += f"  {prefix} return:         {values}\n"

            if self._extra_pdg_ids_to_kill:
                values = ", ".join(
                    _format_pdg_id(pp)
                    for pp in sorted(self._extra_pdg_ids_to_kill)
                )
                final_message += f"  └ don't return:   {values}\n"

            final_message += "\n"

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
        nested_flags += self._return_leaf_flags
        nested_flags += [f"{cc}_cut" for cc in self._cut_definitions]
        nested_flags += [f"include_{pp}" for pp in self._include_processes]
        return nested_flags

    @property
    def _return_leaf_flags(self):
        return [
            f"return_{flag}"
            for flags in self._return_flags.values()
            for flag in flags
        ]

    @property
    def _charged_return_leaf_flags(self):
        return [
            flag
            for flag in self._return_leaf_flags
            if flag.removeprefix("return_")
            not in self._neutral_only_return_flags
        ]

    @property
    def _neutral_only_return_leaf_flags(self):
        return [
            flag
            for flag in self._return_leaf_flags
            if flag.removeprefix("return_")
            in self._neutral_only_return_flags
        ]


    def _get_raw_settings(self):
        veto = self._engine._physics_settings_veto_list
        return {
            name: getattr(self, f"_{name}")
            for name in self._leaf_flags
            if name not in veto
        }

    def _set_raw_settings(self, settings):
        for name, value in settings.items():
            object.__setattr__(self, f"_{name}", value)


    def __getattribute__(self, item):
        # Always use base lookup inside this method
        obj_get = object.__getattribute__
        try:
            engine = obj_get(self, "_engine")
        except AttributeError:
            # _engine does not exist yet during construction
            engine = None

        if item.startswith("return_") and item != "return_pdg_id":
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


    @super_if_being_constructed
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


# ========================
# === Helper functions ===
# ========================

def _format_pdg_id(pdg_id):
    try:
        name = pdg.get_name_from_pdg_id(
            pdg_id,
            long_name=False,
            subscripts=False,
        )
    except ValueError:
        return str(pdg_id)
    else:
        return f"{pdg_id} ({name})"


# ====================================
# === Auto-Generate Getter/Setters ===
# ====================================

def _make_bool_setting(name, default):
    private_name = f"_{name}"
    def getter(self):
        value = getattr(self, private_name)
        if value is None:
            return bool(default(self)) if callable(default) else bool(default)
        return bool(value)
    def setter(self, value):
        if value is not None and not isinstance(value, bool):
            self._error(f"`{name}` has to be a boolean!")
        object.__setattr__(self, private_name, value)
    return property(getter, setter)

def _make_group_property(group, children, prefix=""):
    def getter(self):
        values = [
            getattr(self, f"{prefix}{child}")
            for child in children
        ]
        if all(values):
            return True
        if not any(values):
            return False
        return None
    def setter(self, value):
        if value is not None and not isinstance(value, bool):
            self._error(f"`{prefix}{group}` has to be a boolean!")
        for child in children:
            setattr(self, f"{prefix}{child}", value)
    return property(getter, setter)


for group, flags in PhysicsSettingsHelper._return_flags.items():
    # Return individual flags for each particle type
    for flag, spec in flags.items():
        name = f"return_{flag}"
        default = spec["default"]
        setattr(
            PhysicsSettingsHelper,
            name,
            _make_bool_setting(name, default),
        )

    # Return group flags for each particle group (e.g. return_baryons)
    children = list(flags)
    # For photons and ions, the group itself is already the leaf property.
    if len(children) > 1:
        setattr(
            PhysicsSettingsHelper,
            f"return_{group}",
            _make_group_property(group, children, prefix="return_"),
        )

for modifier, default in PhysicsSettingsHelper._return_modifiers.items():
    name = f"return_{modifier}"
    setattr(
        PhysicsSettingsHelper,
        name,
        _make_bool_setting(name, default),
    )

for modifier, default in PhysicsSettingsHelper._include_processes.items():
    name = f"include_{modifier}"
    setattr(
        PhysicsSettingsHelper,
        name,
        _make_bool_setting(name, default),
    )
