# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

import numpy as np
import scipy.constants as sc

import xobjects as xo

# xo.Strings are complicated
# They come in predefined size of 7, 15, 23, 31, ... characters
# Once the string is assigned for the first time, the max size is defined,
# and extending the string is not possible.
# This means that, e.g. if we create a K2Collimator with Carbon and later change
# that to Aluminium, the name of the latter will be clipped to 'Alumini' (7 chars)

_DEFAULT_NAME = 'NO NAME'
_NAME_MAX_LENGTH = 55  # 55 chars + 1 stop-byte + 8 bytes for length = 64 bytes total

# CLASS INHERITANCE DOES NOT WORK WITH DYNAMIC XOFIELDS
# Because then the memory offsets between parent and child might no longer be the same,
# as xobjects enforces a strict order: static fields first, and then dynamic fields.
# See struct.py, in __new__ of MetaStruct

# A HybridClass needs something to depend on, otherwise the class is added twice in the cdefs during compilation

class Material(xo.HybridClass):
    _xofields = {
        'Z':                        xo.Float64,     # zatom
        'A':                        xo.Float64,     # anuc [g/mol]
        '_ZA_mean':                 xo.Float64,     # [mol of electrons per gram]   electron_density = ZA_mean * density * Avogadro
        'density':                  xo.Float64,     # rho [g/cm^3]
        'radiation_length':         xo.Float64,     # radl [m]
        'excitation_energy':        xo.Float64,     # exenergy [eV]
        '_state':                   xo.Int8,        # 0=solid, 1=liquid, 2=gas, -1=unknown
        '_temperature':             xo.Float64,     # [K]
        '_pressure':                xo.Float64,     # [Pascal]
        # Everest parameters
        'nuclear_radius':           xo.Float64,     # emr
        'nuclear_elastic_slope':    xo.Float64,     # bnref [g/cm^2]
                # slope from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.21.3010
        'cross_section':            xo.Float64[6],  # csref [barn]
                # Index 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
                #       4:Single Diffractive pp or pn, 5:Coulomb for t above mcs
        'hcut':                     xo.Float64,
        '_only_mcs':                xo.Int8,
        # Metadata
        '_name':                    xo.String,
        '_fluka_name':              xo.String,
        '_geant4_name':             xo.String,
    }

    _depends_on = [xo.Float64]
    _skip_in_to_dict  = ['_state', '_temperature', '_pressure', '_name', '_fluka_name', '_geant4_name', '_only_mcs']
    _store_in_to_dict = ['state', 'temperature', 'pressure', 'name', 'fluka_name', 'geant4_name']

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            Z = kwargs.get('Z')
            A = kwargs.get('A')
            if Z is None:
                raise ValueError('Z must be provided')
            if A is None:
                raise ValueError('A must be provided')
            if kwargs.get('density') is None:
                raise ValueError('density must be provided')
            if '_ZA_mean' not in kwargs:
                # Need an if as it can be provided from CompoundMaterial or MixtureMaterial
                kwargs['_ZA_mean'] = Z / A
            if 'radiation_length' not in kwargs:
                density = kwargs['density']
                kwargs['radiation_length'] = _approximate_radiation_length(Z, A, density)
            for kk in ('name', 'geant4_name'):
                to_assign[kk] = kwargs.pop(kk, None)
                kwargs[f'_{kk}'] = 'NO NAME'.ljust(_NAME_MAX_LENGTH)
            to_assign['fluka_name'] = kwargs.pop('fluka_name', None)
            kwargs['_fluka_name'] = 'NO NAME'.ljust(8)
            to_assign['state'] = kwargs.pop('state', None)
            to_assign['temperature'] = kwargs.pop('temperature', None)
            to_assign['pressure'] = kwargs.pop('pressure', None)
            # TODO: this should be better
            kwargs.setdefault('hcut', 0.02)
            kwargs.setdefault('cross_section', [0., 0., 0., 0., 0., 0.])
            kwargs['_only_mcs'] = True
        super().__init__(**kwargs)
        self._frozen = False  # To be able to use the property setters
        for kk, vv in to_assign.items():
            setattr(self, kk, vv)
        if self.full_everest_supported:
            self._only_mcs = False

    def __repr__(self):
        return f"<{self} at {hex(id(self))}>"

    def __str__(self):
        if self.name and (self.name.lower().startswith('liquid') \
        or self.name.lower().startswith('solid') \
        or self.name.lower().startswith('gas')):
            name = self.name
        else:
            state = 'gaseous' if self.state == 'gas' else f'{self.state}'
            name = f"{state} {self.name}" if self.state else self.name
        cls = self.__class__.__name__
        supported_engines = []
        if self.full_everest_supported:
            supported_engines.append('Everest (full)')
        else:
            supported_engines.append('Everest (MCS only)')
        if self.fluka_name is not None:
            supported_engines.append('FLUKA')
        if self.geant4_name is not None:
            supported_engines.append('Geant4')
        if supported_engines:
            supported_engines = ' <' + ', '.join(supported_engines) + '>'
        else:
            supported_engines = ''
        Z = f'{self.Z:.2f}' if not np.isclose(self.Z, int(self.Z), atol=1e-15) else int(self.Z)
        A = f'{self.A:.3f}' if not np.isclose(self.A, int(self.A), atol=1e-15) else f'[{int(self.A)}]'
        return f"{cls}({name}, Z={Z}, A={A}, density={self.density:.4f} g/cm^3)" \
             + f"{supported_engines}"

    def adapt(self, inplace=False, **kwargs):
        thisdict = self.to_dict()
        thisdict.pop('__class__')
        thisdict.pop('fluka_name', None)
        thisdict.pop('geant4_name', None)
        thisdict.update(kwargs)
        if not (np.isclose(self.Z, thisdict['Z']) and np.isclose(self.A, thisdict['A']) \
                and np.isclose(self.density, thisdict['density'])) \
        and kwargs.get('radiation_length') is None:
            thisdict['radiation_length'] = _approximate_radiation_length(thisdict['Z'],
                                                    thisdict['A'], thisdict['density'])
        if inplace:
            for kk, vv in thisdict.items():
                setattr(self, kk, vv)
            return self
        return self.__class__(**thisdict)

    @property
    def electron_density(self):
        return self._ZA_mean * sc.Avogadro # [electrons/g]

    @property
    def plasma_energy(self):
        return np.sqrt(self._ZA_mean * self.density) * 28.816

    @property
    def name(self):
        name = self._name.strip()
        if name == _DEFAULT_NAME:
            return None
        return name

    @name.setter
    def name(self, name):
        self._assert_not_frozen('name')
        from xcoll.materials.database import db as mdb
        old_name = self.name
        if name is None:
            name = _DEFAULT_NAME
            if old_name is not None:
                mdb.remove_material(old_name)
        elif old_name is None:
            mdb[name] = self
        else:
            mdb.rename_material(old_name, name)
        self._name = name

    @property
    def fluka_name(self):
        name = self._fluka_name.strip()
        if name == _DEFAULT_NAME:
            return None
        return name

    @fluka_name.setter
    def fluka_name(self, name):
        self._assert_not_frozen('fluka_name')
        from xcoll.materials.database import db as mdb
        old_name = self.fluka_name
        if name is None:
            name = _DEFAULT_NAME
        else:
            name = name[:8].upper().ljust(8)
            if old_name is None:
                mdb.fluka[name] = self
            else:
                mdb.fluka.rename_material(old_name, name)
        self._fluka_name = name

    @property
    def geant4_name(self):
        name = self._geant4_name.strip()
        if name == _DEFAULT_NAME:
            return None
        return name

    @geant4_name.setter
    def geant4_name(self, name):
        self._assert_not_frozen('geant4_name')
        from xcoll.materials.database import db as mdb
        old_name = self.geant4_name
        if name is None:
            name = _DEFAULT_NAME
        elif old_name is None:
            mdb.geant4[name] = self
        else:
            mdb.geant4.rename_material(old_name, name)
        self._geant4_name = name

    @property
    def full_everest_supported(self):
        return self.excitation_energy != 0 and self.nuclear_radius != 0 \
        and self.nuclear_elastic_slope != 0 and self.hcut != 0 \
        and not all([cc == 0 for cc in self.cross_section])

    @property
    def state(self):
        match self._state:
            case 0:
                return 'solid'
            case 1:
                return 'liquid'
            case 2:
                return 'gas'
            case _:
                return None

    @state.setter
    def state(self, val):
        self._assert_not_frozen('state')
        match val.lower() if isinstance(val, str) else val:
            case 'solid':
                self._state = 0
            case 'liquid':
                self._state = 1
            case 'gas':
                self._state = 2
            case None:
                self._state = -1
            case _:
                raise ValueError("state must be one of 'solid', 'liquid', 'gas', or None")

    @property
    def temperature(self):
        if self._temperature <= 0:
            return None
        return self._temperature

    @temperature.setter
    def temperature(self, val):
        self._assert_not_frozen('temperature')
        if val is None:
            self._temperature = -1
        elif val <= 0:
            raise ValueError('Temperature must be in Kelvin and strictly positive')
        else:
            self._temperature = val

    @property
    def pressure(self):
        if self._pressure <= 0:
            return None
        return self._pressure

    @pressure.setter
    def pressure(self, val):
        self._assert_not_frozen('pressure')
        if val is None:
            self._pressure = -1
        elif val <= 0:
            raise ValueError('Pressure must be in Pascal and strictly positive')
        else:
            self._pressure = val

    def __setattr__(self, name, value):
        if name not in ['_xobject', '_frozen']:
            self._assert_not_frozen(name)
        return super().__setattr__(name, value)

    def _assert_not_frozen(self, name):
        if hasattr(self, '_frozen') and self._frozen:
            raise ValueError(f'Material is frozen, cannot change attribute {name}.')

    # The largest atomic number that can be handled by FLUKA is 100.

    def _gen_geant4_code(self):
        if self.geant4_name:
            raise ValueError(f'Material already has a Geant4 name assigned: {self.geant4_name}.')
        if self.name is None:
            raise ValueError('Material must have a name to generate Geant4 code.')
        name = f'Xcoll_{self.name}'
        code = f"{name} : matdef, Z={self.Z}, A={self.A}, density={self.density}"
        if self.temperature:
            code += f", T={self.temperature}"
        if self.pressure:
            code += f", P={self.pressure}"
        if self.state:
            code += f", state={self.state}"
        code += ";"
        frozen = self._frozen
        self._frozen = False
        self.geant4_name = name
        self._frozen = frozen
        return code


class CrystalMaterial(Material):
    _xofields = Material._xofields | {
        'crystal_nuclear_length':   xo.Float64,     # dlyi
        'crystal_plane_distance':   xo.Float64,     # ai  [mm]
        'crystal_potential':        xo.Float64,     # eum  [eV]
        'nuclear_collision_length': xo.Float64      # collnt [m]
    }

    _depends_on = [Material]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def from_material(cls, material, **kwargs):
        kwargs.setdefault('name', f'{material.name}Crystal')
        thisdict = material.to_dict()
        thisdict.update(kwargs)
        thisdict.pop('__class__')
        thisdict.pop('fluka_name', None)   # Need to define how to deal with these
        thisdict.pop('geant4_name', None)  # Need to define how to deal with these
        return cls(**thisdict)


class CompoundMaterial(Material):
    # Made from atomic elements
    _xofields = Material._xofields | {
        '_components': xo.String[:],
        '_n_atoms':    xo.Int64[:]
    }
    _depends_on = [Material]
    _store_in_to_dict = Material._store_in_to_dict
    _skip_in_to_dict  = Material._skip_in_to_dict

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            components = kwargs.pop('components', None)
            if components is None:
                components = kwargs.pop('_components', None)
                if components is None:
                    raise ValueError('Variable `components` must be provided')
            n_atoms = kwargs.pop('n_atoms', None)
            if n_atoms is None:
                n_atoms = kwargs.pop('_n_atoms', None)
                if n_atoms is None:
                    raise ValueError('Variable `n_atoms` must be provided')
            if len(components) < 2:
                raise ValueError('CompoundMaterial must have at least two components')
            from xcoll.materials.database import db as mdb
            components = [mdb[el] if isinstance(el, str) else el for el in components]
            if any([not isinstance(el, Material) for el in components]):
                raise ValueError('All components must be of type Material')
            if any([isinstance(el, CrystalMaterial) for el in components]):
                raise NotImplementedError("Cannot make CompoundMaterial out of CrystalMaterials")
            if any([isinstance(el, MixtureMaterial) for el in components]):
                raise NotImplementedError("Cannot make CompoundMaterial out of MixtureMaterials")
            if any([nn <= 0 for nn in n_atoms]):
                raise ValueError('All n_atoms must be strictly positive')
            # Resolve nested compounds
            if any([isinstance(el, CompoundMaterial) for el in components]):
                this_components = []
                this_n_atoms = []
                for el, nn in zip(components, n_atoms):
                    if isinstance(el, CompoundMaterial):
                        this_components.extend(el.components)
                        this_n_atoms.extend([n * nn for n in el.n_atoms])
                    else:
                        this_components.append(el)
                        this_n_atoms.append(nn)
                components = this_components
                n_atoms = this_n_atoms
            if any([el.name is None for el in components]):
                raise ValueError('All components must have a name')
            if len(components) != len(n_atoms):
                raise ValueError('Variables `components` and `n_atoms` must have the same length')
            kwargs['_components'] = [el.name.ljust(_NAME_MAX_LENGTH) for el in components]
            kwargs['_n_atoms'] = n_atoms
            if 'Z' in kwargs:
                raise ValueError("Cannot provide Z for CompoundMaterial, it is computed from the components.")
            kwargs['Z'] = sum([el.Z * nn for el, nn in zip(components, n_atoms)]) / sum(n_atoms)
            if 'A' in kwargs:
                raise ValueError("Cannot provide A for CompoundMaterial, it is computed from the components.")
            kwargs['A'] = sum([el.A * nn for el, nn in zip(components, n_atoms)]) / sum(n_atoms)
            if 'density' not in kwargs:
                kwargs['density'] = 0.0
            if 'radiation_length' not in kwargs:
                kwargs['radiation_length'] = 0.0
            if 'excitation_energy' not in kwargs:
                kwargs['excitation_energy'] = 0.0
        super().__init__(**kwargs)
        self._frozen = False
        self._components = [el.strip() for el in self._components]   # Remove padding spaces
        if np.isclose(self.density, 0, atol=1e-15):
            self.density = _inverse_weighted_mean([el.density for el in components],
                                                  self.mass_fractions)
        if np.isclose(self.radiation_length, 0, atol=1e-15):
            self.radiation_length = _inverse_weighted_mean(
                                                [el.radiation_length for el in components],
                                                self.mass_fractions)
        if np.isclose(self.excitation_energy, 0, atol=1e-15) \
        and all([el.excitation_energy > 0 for el in components]):
            fractions  = self.mass_fractions * np.array([el.Z for el in components])
            fractions /= np.array([el.A for el in components])
            self.excitation_energy = _logarithmic_mean(
                                            [el.excitation_energy for el in components],
                                            fractions)
        self._ZA_mean = sum([el.Z * fr / el.A for el, fr in zip(components, self.mass_fractions)])
        if self.full_everest_supported:
            self._only_mcs = False

    @property
    def components(self):
        from xcoll.materials.database import db as mdb
        return [mdb[el] for el in self._components]

    @property
    def n_atoms(self):
        return self._n_atoms

    @property
    def composition(self):
        return [[el.short_name, nn.tolist()] for el, nn in zip(self.components, self.n_atoms)]

    @property
    def molar_mass(self):
        return sum([el.A * nn for el, nn in zip(self.components, self.n_atoms)])

    @property
    def mass_fractions(self):
        return np.array([(el.A * nn) / self.molar_mass for el, nn in zip(self.components, self.n_atoms)])

    def _gen_geant4_code(self):
        if self.geant4_name:
            raise ValueError(f'Material already has a Geant4 name assigned: {self.geant4_name}.')
        if any([el.geant4_name is None for el in self.components]):
            raise ValueError('All components must have a geant4_name to generate Geant4 code.')
        if self.name is None:
            raise ValueError('Material must have a name to generate Geant4 code.')
        name = f'Xcoll_{self.name}'
        code = f"{name} : matdef, density={self.density}"
        if self.temperature:
            code += f", T={self.temperature}"
        if self.pressure:
            code += f", P={self.pressure}"
        if self.state:
            code += f", state={self.state}"
        components = [f'"{el.geant4_name}"' for el in self.components]
        code += f", components=[{','.join(components)}]"
        code += f", componentsWeights={{{','.join([f'{nn}' for nn in self.n_atoms])}}};"
        frozen = self._frozen
        self._frozen = False
        self.geant4_name = name
        self._frozen = frozen
        return code


class MixtureMaterial(Material):
    # Made from atomic elements
    _xofields = Material._xofields | {
        '_components':     xo.String[:],
        '_mass_fractions': xo.Float64[:]
    }
    _depends_on = [Material]
    _store_in_to_dict = Material._store_in_to_dict
    _skip_in_to_dict  = Material._skip_in_to_dict

    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            components = kwargs.pop('components', None)
            if components is None:
                components = kwargs.pop('_components', None)
                if components is None:
                    raise ValueError('Variable `components` must be provided')
            mass_fractions = kwargs.pop('mass_fractions', None)
            if mass_fractions is None:
                mass_fractions = kwargs.pop('_mass_fractions', None)
                if mass_fractions is None:
                    raise ValueError('Variable `mass_fractions` must be provided')
            if len(components) < 2:
                raise ValueError('MixtureMaterial must have at least two components')
            from xcoll.materials.database import db as mdb
            components = [mdb[el] if isinstance(el, str) else el for el in components]
            if any([not isinstance(el, Material) for el in components]):
                raise ValueError('All components must be of type Material')
            if any([isinstance(el, CrystalMaterial) for el in components]):
                raise NotImplementedError("Cannot make MixtureMaterial out of CrystalMaterials")
            if any([nn <= 0 for nn in mass_fractions]):
                raise ValueError('All mass_fractions must be strictly positive')
            # Resolve nested compounds
            if any([isinstance(el, (CompoundMaterial, MixtureMaterial)) for el in components]):
                this_components = []
                this_mass_fractions = []
                for el, fr in zip(components, mass_fractions):
                    if isinstance(el, (CompoundMaterial, MixtureMaterial)):
                        this_components.extend(el.components)
                        this_mass_fractions.extend([ffr * fr for ffr in el.mass_fractions])
                    else:
                        this_components.append(el)
                        this_mass_fractions.append(fr)
                components = this_components
                mass_fractions = this_mass_fractions
            if any([el.name is None for el in components]):
                raise ValueError('All components must have a name')
            if len(components) != len(mass_fractions):
                raise ValueError('Variables `components` and `mass_fractions` must have the same length')
            kwargs['_components'] = [el.name.ljust(_NAME_MAX_LENGTH) for el in components]
            mass_fractions = np.array(mass_fractions) / sum(mass_fractions)
            kwargs['_mass_fractions'] = mass_fractions
            if 'Z' in kwargs:
                raise ValueError("Cannot provide Z for MixtureMaterial, it is computed from the components.")
            kwargs['Z'] = sum([el.Z * fr for el, fr in zip(components, mass_fractions)])
            if 'A' in kwargs:
                raise ValueError("Cannot provide A for MixtureMaterial, it is computed from the components.")
            kwargs['A'] = sum([el.A * fr for el, fr in zip(components, mass_fractions)])
            kwargs['_ZA_mean'] = sum([el.Z * fr / el.A for el, fr in zip(components, mass_fractions)])
            if 'density' not in kwargs:
                kwargs['density'] = _inverse_weighted_mean([el.density for el in components],
                                                           mass_fractions)
            if 'radiation_length' not in kwargs:
                kwargs['radiation_length'] =_inverse_weighted_mean(
                                                [el.radiation_length for el in components],
                                                mass_fractions)
            if 'excitation_energy' not in kwargs:
                fractions  = mass_fractions * np.array([el.Z for el in components])
                fractions /= np.array([el.A for el in components])
                kwargs['excitation_energy'] = _logarithmic_mean(
                                                [el.excitation_energy for el in components],
                                                fractions)
        super().__init__(**kwargs)
        self._frozen = False
        self._components = [el.strip() for el in self._components]   # Remove padding spaces
        if self.full_everest_supported:
            self._only_mcs = False

    @property
    def components(self):
        from xcoll.materials.database import db as mdb
        return [mdb[el] for el in self._components]

    @property
    def mass_fractions(self):
        return self._mass_fractions

    @property
    def composition(self):
        return [[el.short_name, nn.tolist()] for el, nn in zip(self.components, self.mass_fractions)]

    def _gen_geant4_code(self):
        if self.geant4_name:
            raise ValueError(f'Material already has a Geant4 name assigned: {self.geant4_name}.')
        if any([el.geant4_name is None for el in self.components]):
            raise ValueError('All components must have a geant4_name to generate Geant4 code.')
        if self.name is None:
            raise ValueError('Material must have a name to generate Geant4 code.')
        name = f'Xcoll_{self.name}'
        code = f"{name} : matdef, density={self.density}"
        if self.temperature:
            code += f", T={self.temperature}"
        if self.pressure:
            code += f", P={self.pressure}"
        if self.state:
            code += f", state={self.state}"
        components = [f'"{el.geant4_name}"' for el in self.components]
        code += f", components=[{','.join(components)}]"
        code += f", componentsFractions={{{','.join([f'{nn}' for nn in self.mass_fractions])}}};"
        frozen = self._frozen
        self._frozen = False
        self.geant4_name = name
        self._frozen = frozen
        return code

# * Antico
# * ..+....1....+....2....+....3....+....4....+....5....+....6....+....7....+..
# MATERIAL                             2.7                              ANTICO
# COMPOUND      -97.25  ALUMINUM       -.6   SILICON      -0.5      IRONANTICO
# COMPOUND         -.4  MAGNESIU       -.4  MANGANES      -.35  CHROMIUMANTICO
# COMPOUND         -.2  TITANIUM       -.2      ZINC       -.1    COPPERANTICO


class RefMaterial(Material):
    ''' A string to pass a material that exists in Geant4 or FLUKA.'''
    _xofields = Material._xofields
    _depends_on = [Material]
    def __init__(self, **kwargs):
        if '_xobject' not in kwargs:
            if 'name' not in kwargs:
                raise ValueError("ReferenceMaterial must have a name.")
            kwargs = {'_name': kwargs['name']}
        super().__init__(**kwargs)
        self._frozen = True  # ReferenceMaterial is always frozen


def _approximate_radiation_length(Z, A, density):
    # Tsai's more accurate expression (PDG'24 Sec 34.4.2)
    match int(round(Z)):
        case 1:  # Hydrogen
            L_rad = 5.31
            L_rad_prime = 6.144
        case 2:  # Helium
            L_rad = 4.79
            L_rad_prime = 5.621
        case 3:  # Lithium
            L_rad = 4.74
            L_rad_prime = 5.805
        case 4:  # Beryllium
            L_rad = 4.71
            L_rad_prime = 5.924
        case _:
            L_rad = np.log(184.15 * Z ** (-1.0 / 3.0))
            L_rad_prime = np.log(1194.0 * Z ** (-2.0 / 3.0))
    denom = Z**2 * (L_rad - _f_coulomb(Z)) + Z * L_rad_prime
    X0_mass = 716.405 * A / denom  # in g/cm^2
    return X0_mass / density / 100  # in meters

def _f_coulomb(Z):
    # Coulomb correction f(Z) from PDG.
    a = sc.alpha * Z
    a2 = a*a
    return a2 * (
        1/(1+a2)
        + 0.20206
        - 0.0369*a2
        + 0.0083*a2*a2
        - 0.0020*a2*a2*a2
    )

# --- Mixture rules ------------------------------------------------------------

def _inverse_weighted_mean(vals, weights):
    return 1 / (np.array(weights) / np.array(vals)).sum()

def _logarithmic_mean(vals, weights):
    return np.exp(np.sum(np.array(weights) * np.log(np.array(vals))) / np.sum(weights))
