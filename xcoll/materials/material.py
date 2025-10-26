# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import scipy.constants as sc
from collections import defaultdict

import xobjects as xo
from xtrack.line import _dicts_equal

from .parameters import (_approximate_radiation_length, _default_excitation_energies,
                         _combine_radiation_lengths, _combine_excitation_energies,
                         _average_Z_over_A, _effective_Z2, _estimate_density)


_materials_context = xo.ContextCpu()


# CLASS INHERITANCE DOES NOT WORK WITH DYNAMIC XOFIELDS
# Because then the memory offsets between parent and child might no longer be the same,
# as xobjects enforces a strict order: static fields first, and then dynamic fields.
# See struct.py, in __new__ of MetaStruct


class Material(xo.HybridClass):
    _xofields = {
        # Essential fields
        '_density':                 xo.Float64,     # [g/cm^3], can be auto-calculated for compounds
        # Auto-calculated fields
        '_ZA_mean':                 xo.Float64,     # [mol of electrons per gram]
        '_Z2_eff':                  xo.Float64,     # Effective Z for Rutherford scattering
        '_atoms_per_volume':        xo.Float64,     # [atoms/m^3]
        '_num_nucleons_eff':        xo.Float64,     # Effective number of nucleons for nuclear interactions
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
    }

    _depends_on = [xo.Float64] # A HybridClass needs something to depend on, otherwise the class is added twice in the cdefs during compilation

    _skip_in_to_dict  = ['_ZA_mean', '_Z2_eff', '_density', '_radiation_length',
                         '_excitation_energy', '_atoms_per_volume',
                         '_num_nucleons_eff', '_nuclear_radius',
                         '_nuclear_elastic_slope', '_cross_section', '_hcut']
    _store_in_to_dict = ['Z', 'A', 'components', 'n_atoms', 'mass_fractions',
                         'density', 'radiation_length', 'excitation_energy',
                         'nuclear_radius', 'nuclear_elastic_slope',
                         'cross_section', 'hcut', 'state', 'temperature',
                         'pressure', 'info', 'name', 'short_name',
                         'geant4_name', 'fluka_name']


    # ======================
    # === Initialisation ===
    # ======================

    def __init__(self, **kwargs):
        if '_xobject' in kwargs and kwargs['_xobject'] is not None:
            super().__init__(**kwargs)
            return

        # Create xobject with all invalid values (-1)
        xokwargs = kwargs.pop('_xokwargs', {})
        xokwargs['_ZA_mean'] = kwargs.pop('_ZA_mean', -1.)
        xokwargs['_Z2_eff'] = kwargs.pop('_Z2_eff', -1.)
        xokwargs['_radiation_length'] = kwargs.pop('_radiation_length', -1.)
        xokwargs['_excitation_energy'] = kwargs.pop('_excitation_energy', -1.)
        xokwargs['_atoms_per_volume'] = kwargs.pop('_atoms_per_volume', -1.)
        xokwargs['_num_nucleons_eff'] = kwargs.pop('_num_nucleons_eff', -1.)
        xokwargs['_density'] = kwargs.get('_density', -1.)
        xokwargs['_nuclear_radius'] = kwargs.pop('_nuclear_radius', -1.)
        xokwargs['_nuclear_elastic_slope'] = kwargs.pop('_nuclear_elastic_slope', -1.)
        xokwargs['_cross_section'] = kwargs.pop('_cross_section', [-1., -1., -1., -1., -1., -1.])
        xokwargs['_hcut'] = kwargs.pop('_hcut', -1.)
        xokwargs['_context'] = kwargs.pop('_context', _materials_context)  # This is needed to get all materials in the same buffer (otherwise Xtrack tests fail)
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

        # Calculate dependent properties
        self.density = kwargs.pop('density', None)                     # Can be provided for more precision (obligatory for elements)
        self.radiation_length = kwargs.pop('radiation_length', None)   # Can be provided for more precision
        self.excitation_energy = kwargs.pop('excitation_energy', None) # Can be provided for more precision
        self.update_vars()

        # Assign optional properties
        for kk in ['nuclear_radius', 'nuclear_elastic_slope', 'cross_section', 'hcut']:
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
        if 'density' not in kwargs:
            raise ValueError('density must be provided for an elemental Material')


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
        if any([not isinstance(el, Material) for el in components]):
            raise ValueError('All components must be of type Material')
        if any([isinstance(el, CrystalMaterial) for el in components]):
            raise NotImplementedError("Cannot make Material out of CrystalMaterials")
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
            raise ValueError("One of `n_atoms`, `mass_fractions` or "
                             "`volume_fractions` must be provided")

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
        mass_fractions = np.array([mf * el.molar_mass for el, mf in
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
            elif this_cls == 'CrystalMaterial':
                return CrystalMaterial.from_dict(dct)
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
            if mdb[name] != mat:
                print(f"Warning: Material {name} is out of sync with database.")
                mat._out_of_sync = True
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
        return dct

    def adapt(self, inplace=False, **kwargs):
        thisdict = self.to_dict()
        thisdict.pop('__class__')
        thisdict.pop('name', None)
        thisdict.pop('short_name', None)
        thisdict.pop('fluka_name', None)
        thisdict.pop('geant4_name', None)
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
        if self.full_everest_supported:
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
        return _dicts_equal(self._xobject._to_dict(),
                            other._xobject._to_dict())

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
        if self.components is not None:
            raise ValueError('Cannot set Z for a compound material')
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
        if self.components is not None:
            raise ValueError('Cannot set A for a compound material')
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
        else:
            # Use molar_mass such that this can be used for nested compounds
            mf = np.array([fr / el.molar_mass for el, fr
                                in zip(self.components, self.mass_fractions)])
            return mf / mf.sum()

    @property
    def atomic_fractions(self):
        return self.molar_fractions

    @property
    def volume_fractions(self):
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
            if self.components is not None:
                self._density = _estimate_density(self.components,
                                                  self.mass_fractions)
            else:
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
            if self.components is not None:
                self._radiation_length = _combine_radiation_lengths(
                            self.components, self.mass_fractions, self.density)
            else:
                if self.Z is not None and self.A is not None \
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
    def electron_density(self):
        return self._ZA_mean * sc.Avogadro # [electrons/g]

    @property
    def plasma_energy(self):
        return np.sqrt(self._ZA_mean * self.density) * 28.816

    @property
    def atoms_per_volume(self):
        # TODO: do we want to use a different A? Maybe an average instead of molar mass?
        mass = self.molar_mass or self.average_molar_mass
        if self.density and mass:
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
        if A1_3:
            return 2*np.pi*dn*(3/4/np.pi)**(1/3) * A1_3


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

# * Antico
# * ..+....1....+....2....+....3....+....4....+....5....+....6....+....7....+..
# MATERIAL                             2.7                              ANTICO
# COMPOUND      -97.25  ALUMINUM       -.6   SILICON      -0.5      IRONANTICO
# COMPOUND         -.4  MAGNESIU       -.4  MANGANES      -.35  CHROMIUMANTICO
# COMPOUND         -.2  TITANIUM       -.2      ZINC       -.1    COPPERANTICO

    def _gen_geant4_code(self):
        if self.geant4_name:
            raise ValueError(f'Material already has a Geant4 name assigned: {self.geant4_name}.')
        if self.name is None:
            raise ValueError('Material must have a name to generate Geant4 code.')
        if self.components is not None:
            if any([el.geant4_name is None for el in self.components]):
                raise ValueError('All components must have a geant4_name to generate Geant4 code.')
        name = f'Xcoll_{self.name}'
        if self.components is not None:
            code = f"{name} : matdef, density={self.density}"
        else:
            code = f"{name} : matdef, Z={self.Z}, A={self.A}, density={self.density}"
        if self.temperature:
            code += f", T={self.temperature}"
        if self.pressure:
            code += f", P={self.pressure}"
        if self.state:
            code += f", state={self.state}"
        if self.components is not None:
            components = [f'"{el.geant4_name}"' for el in self.components]
            code += f", components=[{','.join(components)}]"
            if self.n_atoms is not None:
                code += f", componentsWeights={{{','.join([f'{nn}' for nn in self.n_atoms])}}};"
            else:
                code += f", componentsFractions={{{','.join([f'{nn}' for nn in self.mass_fractions])}}}"
        code += ";"
        frozen = self._frozen
        self._frozen = False
        self.geant4_name = name
        self._frozen = frozen
        return code


class CrystalMaterial(Material):
    _xofields = Material._xofields | {
        '_crystal_plane_distance':   xo.Float64,     # ai  [mm]
        '_crystal_potential':        xo.Float64,     # eum  [eV]
        '_nuclear_collision_length': xo.Float64,     # collnt [m]
        '_eta':                      xo.Float64
    }

    _depends_on = [Material]
    _skip_in_to_dict  = Material._skip_in_to_dict + [
                         '_crystal_plane_distance', '_crystal_potential',
                         '_nuclear_collision_length', '_eta']
    _store_in_to_dict = Material._store_in_to_dict + [
                         'crystal_plane_distance', 'crystal_potential',
                         'nuclear_collision_length', 'eta']

    def __init__(self, **kwargs):
        if '_xobject' in kwargs and kwargs['_xobject'] is not None:
            super().__init__(**kwargs)
            return
        xokwargs = {}
        for kk in ('crystal_plane_distance',
                   'crystal_potential', 'eta',
                   'nuclear_collision_length'):
            xokwargs[f'_{kk}'] = kwargs.pop(kk, -1.)
        kwargs['_xokwargs'] = xokwargs
        super().__init__(**kwargs)

    @classmethod
    def from_material(cls, material, **kwargs):
        kwargs.setdefault('name', f'{material.name}Crystal')
        thisdict = material.to_dict()
        thisdict.update(kwargs)
        thisdict.pop('__class__')
        thisdict.pop('short_name', None)   # Need to define how to deal with these
        thisdict.pop('fluka_name', None)   # Need to define how to deal with these
        thisdict.pop('geant4_name', None)  # Need to define how to deal with these
        return cls(**thisdict)

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



class RefMaterial(Material):
    ''' A string to pass a material that exists in Geant4 or FLUKA.'''
    _xofields = Material._xofields
    _depends_on = [Material]
    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            if 'name' not in kwargs:
                raise ValueError("ReferenceMaterial must have a name.")
            name = kwargs.pop('name')
            # Fake values to pass initialisation
            kwargs['Z'] = 1
            kwargs['A'] = 1
            kwargs['density'] = 1
        super().__init__(**kwargs)
        # Invalidate all properties
        self.invalidate()
        # Set name only now to avoid syncing with database
        self._name = name
        self._frozen = True  # ReferenceMaterial is always frozen

    def __str__(self):
        return f"RefMaterial '{self.name}'"


_DEFAULT_MATERIAL = Material(Z=1, A=1, density=1)
_DEFAULT_MATERIAL.invalidate()
_DEFAULT_CRYSTALMATERIAL = CrystalMaterial(Z=1, A=1, density=1)
_DEFAULT_CRYSTALMATERIAL.invalidate()
