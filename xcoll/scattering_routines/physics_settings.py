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
    _global_return_flags = ['all', 'all_charged', 'none']
    _return_flags = {
        'photons': {
            'photons': lambda self: self.ref_is_photon,
        },
        'leptons': {
            'electrons': lambda self: self.ref_is_lepton,
            'muons': lambda self: self.ref_is_lepton,
            'tauons': lambda self: self.ref_is_lepton,
            'neutrinos': lambda self: self.ref_is_lepton,
        },
        'baryons': {
            'protons': lambda self: self.ref_is_baryon or self.ref_is_ion,
            'neutrons': lambda self: self.ref_is_baryon or self.ref_is_ion,
            'other_baryons': lambda self: self.ref_is_baryon,
        },
        'mesons': {
            'pions': lambda self: self.ref_is_meson,
            'kaons': lambda self: self.ref_is_meson,
            'other_mesons': lambda self: self.ref_is_meson,
        },
        'ions': {
            'ions': lambda self: self.ref_is_ion,
        },
    }
    _return_modifiers = {
        'neutral': lambda self: self.ref_is_neutral,
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
        'showers': True,
        'single_coulomb': True,
        'multiple_coulomb': True,
        'ionisation_fluctuations': True,
        'pair_production': True,
        'bremsstrahlung': True,
        'elastic': True,
        'inelastic': True
    }


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
        nested_flags += self._return_leaf_flags()
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
    def ref_is_proton(self):
        return pdg.is_proton(self.ref_id)

    @property
    def ref_is_ion(self):
        return pdg.is_ion(self.ref_id)

    @property
    def ref_is_lepton(self):
        return pdg.is_lepton(self.ref_id)

    @property
    def ref_is_meson(self):
        # TODO: Use pdg.is_meson() once it is implemented in xtrack
        return bool(_is_meson(self.ref_id))

    @property
    def ref_is_baryon(self):
        # TODO: Use pdg.is_baryon() once it is implemented in xtrack
        return bool(_is_baryon(self.ref_id))

    @property
    def ref_is_photon(self):
        # TODO: Use pdg.is_photon() once it is implemented in xtrack
        return abs(self.ref_id) == 22

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


    # =======================
    # === Private Methods ===
    # =======================

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


# ========================
# === Helper functions ===
# ========================

# TODO: Use pdg.is_meson() once it is implemented in xtrack
def _is_meson(pdg_id):
    pid = np.abs(np.asarray(pdg_id))
    q1 = (pid // 100) % 10
    q2 = (pid // 10) % 10
    q3 = (pid // 1000) % 10
    return (
        (pid < 1_000_000_000)
        & (q3 == 0)
        & (q1 >= 1) & (q1 <= 6)
        & (q2 >= 1) & (q2 <= 6)
    )

# TODO: Use pdg.is_baryon() once it is implemented in xtrack
def _is_baryon(pdg_id):
    pid = np.abs(np.asarray(pdg_id))
    q1 = (pid // 1000) % 10
    q2 = (pid // 100) % 10
    q3 = (pid // 10) % 10
    return (
        (pid < 1_000_000_000)
        & (q1 >= 1) & (q1 <= 6)
        & (q2 >= 1) & (q2 <= 6)
        & (q3 >= 1) & (q3 <= 6)
    )


# ====================================
# === Auto-Generate Getter/Setters ===
# ====================================

def _make_bool_setting(name, default):
    private_name = f"_{name}"
    def getter(self):
        value = getattr(self, private_name)
        if value is None:
            return default(self) if callable(default) else default
        return value
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
    for flag, default in flags.items():
        name = f"return_{flag}"
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
