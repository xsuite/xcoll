# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import time
import numpy as np

try:
    from xaux import FsPath, ranID  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath, ranID

from .environment import format_fluka_float
from ... import json
from ...beam_elements.base import BaseCollimator
from ...materials import Material
from ...compare import deep_equal


class FlukaPrototype:
    # This is a registry to keep track of prototypes. If one is already defined,
    # we do not create a new instance but return the existing one
    _registry = []

    def __new__(cls, fedb_series=None, fedb_tag=None, **kwargs):
        _is_null = False
        # If the prototype is already defined, return the existing instance
        if fedb_series is None and fedb_tag is None:
            _is_null = True
        elif fedb_series is None or fedb_tag is None:
            raise ValueError("Both 'fedb_series' and 'fedb_tag' must be provided.")
        else:
            for prototype in FlukaPrototype._registry:
                if prototype.fedb_series.upper() == fedb_series.upper() \
                and prototype.fedb_tag.upper() == fedb_tag.upper() \
                and prototype.__class__ == cls:
                    return prototype
        # Register the new prototype
        self = object.__new__(cls)
        self._is_null = _is_null
        if not _is_null:
            FlukaPrototype._registry.append(self)
        return self

    def __init__(self, fedb_series=None, fedb_tag=None, *, angle=0, side=None, width=None,
                 height=None, length=None, material=None, info=None, extra_commands=None,
                 is_crystal=False, bending_radius=None, _allow_generic=False, is_broken=False,
                 _force_init=False, **kwargs):
        if getattr(self, "_initialized", False) and not _force_init:
            return
        self._idx = None
        self._type = self.__class__.__name__[5:].lower()
        self._file_is_valid = None
        # while True:
        #     self._hash = ranID(length=256)
        #     if self._hash not in [pp._hash for pp in FlukaPrototype._registry if pp is not self]:
        #         break
        if self._is_null:
            self._fedb_series = None
            self._fedb_tag = None
            self._name = None
            self._side = None
            self._angle = None
            self._length = None
            self._width = None
            self._height = None
            self._material = None
            self._is_crystal = None
            self._bending_radius = None
            self._info = None
            self._extra_commands = None
            self._is_broken = None
            self._initialized = True
            return
        if fedb_series == 'generic' and not _allow_generic:
            this_type = self.__class__.__name__[5:].lower()
            raise ValueError("Cannot use 'generic' as fedb_series, unless creating a generic " \
                          + f"{this_type}. Please use xcoll.fluka.create_generic_{this_type}() " \
                          + f"instead.")
        self._fedb_series = fedb_series
        self._fedb_tag = fedb_tag
        self._name = fedb_tag
        if side is not None:
            BaseCollimator.side.fset(self, side)  # This will overwrite the side in the FlukaCollimator
        else:
            self._side = None
        self._angle = angle        # This will adapt the transformation used in the LineBuilder
        self._length = length      # This will adapt the Drift used in pre-tracking (length_front, length_back)
        self._width = width
        self._height = height
        if isinstance(material, dict):
            if material.get('fluka_name') is not None \
            and material['fluka_name'].startswith('XCOLL'):
                # Remove fluka_name for previously auto-generated material,
                # such that the material code can be re-generated
                material.pop('fluka_name')
            material = Material.from_dict(material)
        self._material = material
        self._is_crystal = is_crystal
        self._bending_radius = bending_radius
        self._info = info
        self._extra_commands = extra_commands
        self._is_broken = is_broken
        self._generic_body_file = None
        self._generic_region_file = None
        self._generic_material_file = None
        self._initialized = True

    def __repr__(self):
        if self._is_null:
            return ''
        info = f" ({self.info})" if self.info else ''
        cry = "crystal " if self.is_crystal else ''
        defunct = " <defunct>" if self.is_defunct() else ''
        return f"{self.__class__.__name__} {cry}'{self.name}': " \
             + f"tag {self.fedb_tag} in {self.fedb_series} series{info}{defunct}"

    def __str__(self):
        if self._is_null:
            return ''
        return self.__repr__()

    # def __eq__(self, other):
    #     if not isinstance(other, FlukaPrototype):
    #         return False
    #     if self._is_null and other._is_null:
    #         return True
    #     if self._is_null != other._is_null:
    #         return False
    #     return deep_equal(self.to_dict(), other.to_dict())

    # def __hash__(self):
    #     return self._hash

    def is_defunct(self):
        if self._is_null:
            return False
        # Needed for proper handling of prototypes that are deleted
        return self not in self._registry

    def is_generic(self):
        if self._is_null:
            return False
        return self.fedb_series.lower() == 'generic'

    def delete(self, **kwargs):
        import xcoll as xc
        if self._is_null:
            return
        while self in self._registry:
            self._registry.remove(self)
        # Remove all files associated with the prototype
        for file in self.files:
            if file.exists() or file.is_symlink():
                try:
                    file.unlink()
                except FileNotFoundError:
                    pass
        fedb = xc.fluka.environment.fedb
        meta = fedb / "metadata" / f'{self.fedb_series}_{self.fedb_tag}.bodies.json'
        if meta.exists() or meta.is_symlink():
            try:
                meta.unlink()
            except FileNotFoundError:
                pass

    def to_dict(self):
        if self._is_null:
            return {'__class__': self.__class__.__name__}
        if self.is_defunct():
            raise ValueError(f"Cannot serialize defunct {self._type} '{self.name}'!")
        return {
            '__class__': self.__class__.__name__,
            'name': self.name,
            'fedb_series': self.fedb_series,
            'fedb_tag': self.fedb_tag,
            'side': self.side,
            'angle': self.angle,
            'length': self.length,
            'width': self.width,
            'height': self.height,
            'material': self.material.to_dict() if hasattr(self.material, 'to_dict') else self.material,
            'is_crystal': self.is_crystal,
            'bending_radius': self.bending_radius,
            'info': self.info,
            'extra_commands': self.extra_commands,
            'is_broken': self.is_broken,
        }

    @classmethod
    def from_dict(cls, data):
        cls = data.pop('__class__', None)
        if cls == 'FlukaPrototype':
            return FlukaPrototype(_allow_generic=True, _force_init=True, **data)
        elif cls == 'FlukaAssembly':
            return FlukaAssembly(_allow_generic=True, _force_init=True, **data)
        else:
            raise ValueError(f"Invalid data format for {cls}.")

    @classmethod
    def from_json(cls, path):
        return cls.from_dict(json.load(path))

    @property
    def name(self):
        if self._is_null:
            return None
        return self._name

    @name.setter
    def name(self, val):
        if self._is_null:
            return
        if self.is_defunct():
            raise ValueError(f"Cannot rename defunct {self._type} '{self.name}'!")
        if len(val) > 8:
            raise ValueError(f"Prototype name '{val}' is too long (max. 8 characters).")
        self._name = val

    @property
    def fedb_series(self):
        if self._is_null:
            return None
        return self._fedb_series

    @property
    def fedb_tag(self):
        if self._is_null:
            return None
        return self._fedb_tag

    @property
    def body_file(self):
        import xcoll as xc
        if self._is_null:
            return None
        if self.is_generic():
            return self._generic_body_file
        fedb = xc.fluka.environment.fedb
        file = fedb / "bodies" / f"{self.fedb_series}_{self.fedb_tag}.bodies"
        return file.resolve()

    @body_file.setter
    def body_file(self, path):
        import xcoll as xc
        if self._is_null:
            raise ValueError("Cannot set body_file for a null prototype!")
        if self.is_defunct():
            raise ValueError(f"Cannot set body_file for defunct {self._type} '{self.name}'!")
        if self.is_generic():
            raise ValueError(f"Cannot set body_file for generic {self._type} '{self.name}'! "
                            + "These are generated automatically.")
        path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        fedb = xc.fluka.environment.fedb
        target = fedb / "bodies" / f"{self.fedb_series}_{self.fedb_tag}.bodies"
        if path != target:
            path.copy_to(target, method='mount')
        file = fedb / "metadata" / f'{self.fedb_series}_{self.fedb_tag}.bodies.json'
        json.dump(self.to_dict(), file, indent=2)

    @property
    def material_file(self):
        import xcoll as xc
        if self._is_null:
            return None
        if self.is_generic():
            return self._generic_material_file
        fedb = xc.fluka.environment.fedb
        file = fedb / "materials" / f"{self.fedb_series}_{self.fedb_tag}.assignmat"
        return file.resolve()

    @material_file.setter
    def material_file(self, path):
        import xcoll as xc
        if self._is_null:
            raise ValueError("Cannot set material_file for a null prototype!")
        if self.is_defunct():
            raise ValueError(f"Cannot set material_file for defunct {self._type} '{self.name}'!")
        if self.is_generic():
            raise ValueError(f"Cannot set material_file for generic {self._type} '{self.name}'! "
                            + "These are generated automatically.")
        path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        fedb = xc.fluka.environment.fedb
        path.copy_to(fedb / "materials" / f"{self.fedb_series}_{self.fedb_tag}.assignmat",
                     method='mount')

    @property
    def region_file(self):
        import xcoll as xc
        if self._is_null:
            return None
        if self.is_generic():
            return self._generic_region_file
        fedb = xc.fluka.environment.fedb
        file = fedb / "regions" / f"{self.fedb_series}_{self.fedb_tag}.regions"
        return file.resolve()

    @region_file.setter
    def region_file(self, path):
        import xcoll as xc
        if self._is_null:
            raise ValueError("Cannot set region_file for a null prototype!")
        if self.is_defunct():
            raise ValueError(f"Cannot set region_file for defunct {self._type} '{self.name}'!")
        if self.is_generic():
            raise ValueError(f"Cannot set region_file for generic {self._type} '{self.name}'! "
                            + "These are generated automatically.")
        path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        fedb = xc.fluka.environment.fedb
        path.copy_to(fedb / "regions" / f"{self.fedb_series}_{self.fedb_tag}.regions",
                     method='mount')

    def populate_into_temp_fedb(self, fedb):
        fedb = FsPath(fedb).resolve()
        if not fedb.exists():
            fedb.mkdir(parents=True)
            (fedb / 'assemblies').mkdir(parents=True)
            (fedb / 'bodies').mkdir(parents=True)
            (fedb / 'regions').mkdir(parents=True)
            (fedb / 'materials').mkdir(parents=True)
            (fedb / 'stepsizes').mkdir(parents=True)
        if self.is_generic():
            for f in self.files:
                if f is None or not f.exists():
                    raise ValueError(f"Generic prototype '{self.name}' is missing files!")
            return  # Generic prototypes are generated automatically
        link = fedb / 'bodies' / self.body_file.name
        if not link.exists():
            link.symlink_to(self.body_file)
        mat_link = fedb / 'materials' / f'{self.body_file.stem}.assignmat'
        if not mat_link.exists() and self.material_file is not None:
            mat_link.symlink_to(self.material_file)
        reg_link = fedb / 'regions' / f'{self.body_file.stem}.regions'
        if not reg_link.exists() and self.region_file is not None:
            reg_link.symlink_to(self.region_file)

    @property
    def files(self):
        return [self.body_file, self.material_file, self.region_file]

    def exists(self):
        if self.is_generic():
            return True
        if self._is_null or not self.files or self.is_defunct():
            return False
        return np.all([ff.exists() for ff in self.files])

    def assert_exists(self):
        if not self.exists():
            raise ValueError(f"{self._type.capitalize()} '{self.name}' "
                            + f"does not exist in the FEDB!")

    @property
    def side(self):
        if self._is_null:
            return None
        return BaseCollimator.side.fget(self)

    @property
    def angle(self):
        if self._is_null:
            return None
        return self._angle

    @property
    def length(self):
        if self._is_null:
            return None
        return self._length

    @property
    def width(self):
        if self._is_null:
            return None
        return self._width

    @property
    def height(self):
        if self._is_null:
            return None
        return self._height

    @property
    def material(self):
        if self._is_null:
            return None
        return self._material

    @property
    def is_crystal(self):
        if self._is_null:
            return None
        return self._is_crystal

    @property
    def bending_radius(self):
        if self._is_null:
            return None
        return self._bending_radius

    @property
    def info(self):
        if self._is_null:
            return None
        return self._info

    @property
    def extra_commands(self):
        if self._is_null:
            return None
        return self._extra_commands

    @property
    def is_broken(self):
        if self._is_null:
            return None
        return self._is_broken

    @property
    def fluka_position(self):
        if self._idx is None:
            return None
        # Maximum positions for the parking region:
        # x in [-3000.0, 3000.0], y in [-4000.0, -2000.0], z in [0.0, 1.E5]
        return 0., 0., 0., (self._idx%5-2)*500.0 , -3000., (self._idx//5+1)*1000.0

    def generate_code(self, idx=0, elements=[]):
        import xcoll as xc
        if self._is_null:
            raise ValueError(f"Cannot generate code for null {self._type}!")
        if self.is_defunct():
            raise ValueError(f"Cannot generate code for defunct {self._type} '{self.name}'!")
        if not self.exists():
            raise ValueError(f"{self._type.capitalize()} '{self.name}' "
                           + f"does not exist in the FEDB!")
        self.check_file_valid()
        _type = 'ASSEMBLY' if isinstance(self, FlukaAssembly) else 'PROTOTYPE'
        prot  = f"{_type:9}     {self.name}\n"
        prot += f"FEDB_SERIES   {self.fedb_series}\n"
        prot += f"FEDB_TAG      {self.fedb_tag}\n"
        prot += f"ROT-DEFI  "
        self._idx = idx  # Store the index for fluka_position property
        for value in self.fluka_position:
            prot += format_fluka_float(value)
        prot += " proto\n"
        if self._extra_commands:
            if hasattr(self._extra_commands, "__iter__") \
            and not isinstance(self._extra_commands, str):
                prot += "\n".join(self._extra_commands) + "\n"
            else:
                prot += self._extra_commands + "\n"
        maps = []
        for element in elements:
            if not hasattr(element, 'name'):
                raise ValueError(f"Element {element} has no `name` variable! "
                               + f"Cannot assign to a {self._type}!")
            if not isinstance(element, xc.fluka.engine._element_classes):
                raise ValueError(f"Element {element.name} is not a FLUKA element! "
                               + f"Cannot assign to a {self._type}!")
            if not element.assembly is self:
                raise ValueError(f"Element '{element.name}' is not assigned to "
                               + f"{self._type} '{self.name}'!")
            maps.append(f" {element.name.upper()}")
        if len(maps) > 0:
            prot_map = [f"MAP_ENTRIES  "]
            for element in maps:
                if len(prot_map[-1] + element) > 75:
                    prot_map.append(f"MAP_ENTRIES  ")
                prot_map[-1] += element
            prot += "\n".join(prot_map) + "\n"
        prot += '#'
        return prot

    def _check_name_clash(self, existing_names):
        if self.name.upper() in existing_names:
            i = 0
            while True:
                new_name = f"{self.name}{i}"
                if new_name.upper() not in existing_names:
                    print(f"Warning: {self._type.capitalize()} name {self.name} "
                        + f"already in use. Renaming to '{new_name}'.")
                    self.name = new_name
                    break
                i += 1
        existing_names.append(self.name.upper())

    @classmethod
    def _group_elements_by_prototype(cls, elements):
        import xcoll as xc
        # Sort elements by their assigned prototype/assembly
        all_prototypes = {}
        all_assemblies = {}
        if not hasattr(elements, '__iter__'):
            elements = [elements]
        for element in elements:
            if not isinstance(element, xc.fluka.engine._element_classes):
                raise ValueError(f"Element {element.name} is not a FLUKA element! "
                               + f"Cannot generate prototypes file!")
            assm = element.assembly
            if assm is None:
                raise ValueError(f"Element '{element.name}' has no assigned prototype! "
                               + f"Cannot generate prototypes file!")
            if isinstance(assm, FlukaAssembly):
                if assm not in all_assemblies:
                    all_assemblies[assm] = []
                all_assemblies[assm].append(element)
            elif isinstance(assm, FlukaPrototype):
                if assm not in all_prototypes:
                    all_prototypes[assm] = []
                all_prototypes[assm].append(element)
            else:
                raise ValueError(f"Element '{element.name}' has an invalid assigned prototype! "
                               + f"Cannot generate prototypes file!")
        return all_prototypes, all_assemblies

    @classmethod
    def make_prototypes_file(cls, elements, save=True, path=None):
        all_prototypes, all_assemblies = cls._group_elements_by_prototype(elements)
        # Generate code
        file_content = ["#...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8"]
        i = 0  # Initialise i in case there are no prototypes
        existing_names = []
        # First the prototypes
        for i, (pro, els) in enumerate(all_prototypes.items()):
            assert isinstance(pro, FlukaPrototype)
            pro._check_name_clash(existing_names)
            file_content.append(pro.generate_code(idx=i, elements=els))
        # Then the assemblies
        for j, (assm, els) in enumerate(all_assemblies.items()):
            assert isinstance(assm, FlukaAssembly)
            assm._check_name_clash(existing_names)
            file_content.append(assm.generate_code(idx=i+j+1, elements=els))
        # Sanity check
        all_names = [assm.name for assm in (all_prototypes | all_assemblies).keys()]
        if len(all_names) != len(set(all_names)):
            raise ValueError("Duplicate prototype/assembly names found! "
                           + "Cannot generate prototypes file!")
        # Finalise
        file_content = "\n".join(file_content)
        if save:
            if path is None:
                path = FsPath.cwd()
            with (path / "prototypes.lbp").open("w") as fp:
                fp.write(file_content)
        return file_content

    @classmethod
    def inspect_prototypes_file(cls, elements, prototypes_file=None):
        if prototypes_file is None:
            prototypes_file = FsPath.cwd() / "prototypes.lbp"
        prototypes_in_file = {}
        with FsPath(prototypes_file).open("r") as fp:
            name = None
            for line in fp.readlines():
                line = line.strip()
                if line.startswith("#") or len(line) == 0:
                    continue
                if line.startswith("PROTOTYPE"):
                    name = line.split()[1]
                    prototypes_in_file[name] = {'type': 'PROTOTYPE'}
                elif line.startswith("ASSEMBLY"):
                    name = line.split()[1]
                    prototypes_in_file[name] = {'type': 'ASSEMBLY'}
                elif line.startswith("FEDB_SERIES"):
                    if name is None:
                        raise ValueError(f"File {prototypes_file} malformed ("
                                    + "FEDB_SERIES without PROTOTYPE or ASSEMBLY)!")
                    prototypes_in_file[name]['fedb_series'] = line.split()[1]
                elif line.startswith("FEDB_TAG"):
                    if name is None:
                        raise ValueError(f"File {prototypes_file} malformed ("
                                    + "FEDB_TAG without PROTOTYPE or ASSEMBLY)!")
                    prototypes_in_file[name]['fedb_tag'] = line.split()[1]
                elif line.startswith("MAP_ENTRIES"):
                    if name is None:
                        raise ValueError(f"File {prototypes_file} malformed ("
                                    + "MAP_ENTRIES without PROTOTYPE or ASSEMBLY)!")
                    if 'elements' not in prototypes_in_file[name]:
                        prototypes_in_file[name]['elements'] = []
                    prototypes_in_file[name]['elements'] += line.split()[1:]
        prototypes, assemblies = cls._group_elements_by_prototype(elements)
        for prototype, els in (prototypes | assemblies).items():
            if prototype.name not in prototypes_in_file:
                raise ValueError(f"Prototype {prototype.name} (in {prototype.fedb_series}) not found in prototypes file!")
            else:
                _type = 'ASSEMBLY' if isinstance(prototype, FlukaAssembly) else 'PROTOTYPE'
                if _type != prototypes_in_file[prototype.name]['type']:
                    raise ValueError(f"Wrong type for {prototype.name} (in {prototype.fedb_series}) (expected "
                                  + f"{prototype._type.upper()}, got "
                                  + f"{prototypes_in_file[prototype.name]['type']})!")
                if 'fedb_series' not in prototypes_in_file[prototype.name]:
                    raise ValueError(f"FEDB_SERIES for {prototype.name} (in {prototype.fedb_series}) not found in prototypes file!")
                if prototype.fedb_series != prototypes_in_file[prototype.name]['fedb_series']:
                    raise ValueError(f"Wrong fedb_series for {prototype.name} (in {prototype.fedb_series}) (expected "
                                   + f"{prototype.fedb_series}, got "
                                   + f"{prototypes_in_file[prototype.name]['fedb_series']})!")
                if 'fedb_tag' not in prototypes_in_file[prototype.name]:
                    raise ValueError(f"FEDB_TAG for {prototype.name} (in {prototype.fedb_series}) not found in prototypes file!")
                if prototype.fedb_tag != prototypes_in_file[prototype.name]['fedb_tag']:
                    raise ValueError(f"Wrong fedb_tag for {prototype.name} (in {prototype.fedb_series}) (expected "
                                   + f"{prototype.fedb_tag}, got "
                                   + f"{prototypes_in_file[prototype.name]['fedb_tag']})!")
                if 'elements' not in prototypes_in_file[prototype.name]:
                    raise ValueError(f"MAP_ENTRIES for {prototype.name} (in {prototype.fedb_series}) not found in prototypes file!")
                for element in els:
                    if element.name.upper() not in prototypes_in_file[prototype.name]['elements']:
                        raise ValueError(f"Element {element.name} not found in prototypes file!")

    @property
    def dependant_assemblies(self):
        """Return a list of all assemblies that depend on this prototype."""
        if self._is_null:
            return []
        assemblies = []
        for assembly in FlukaPrototype._registry:
            if isinstance(assembly, FlukaAssembly) \
            and assembly.fedb_series == self.fedb_series:
                if self in assembly.prototypes:
                    assemblies.append(assembly)
        return assemblies

    def check_file_valid(self, raise_error=True):
        return True

    def view(self, show=True, keep_files=False):
        import xcoll as xc
        if self.exists():
            xc.fluka.environment.test_assembly(self.fedb_series, self.fedb_tag, show=show,
                                               keep_files=keep_files)


class FlukaAssembly(FlukaPrototype):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._generic_assembly_file = None

    def delete(self, _ignore_files=False, **kwargs):
        import xcoll as xc
        if self._is_null:
            return
        # Remove prototypes if no other assembly depends on them
        if self.assembly_file.exists() or not _ignore_files:
            try:
                to_delete = []  # Need to do like this to avoid prototypes spawning again (during checking of dependant assemblies)
                for pro in self.prototypes:
                    if pro.dependant_assemblies == [self]:
                        to_delete.append(pro)
                for pro in to_delete:
                    # Only delete the prototype if no other assembly depends on it
                    try:
                        pro.delete()
                    except Exception as e:
                        print(f"Could not remove prototype '{pro.name}': {e}")
            except Exception as e:
                print(f"Could not remove dependent prototypes: {e}")
        # Remove the assembly from the registry of all prototypes
        while self in FlukaPrototype._registry:
            FlukaPrototype._registry.remove(self)
        # Remove all files associated with the assembly
        if self.assembly_file.exists() or self.assembly_file.is_symlink():
            try:
                self.assembly_file.unlink()
            except FileNotFoundError:
                pass
        fedb = xc.fluka.environment.fedb
        meta = fedb / "metadata" / f'{self.fedb_series}_{self.fedb_tag}.lbp.json'
        if meta.exists() or meta.is_symlink():
            meta.unlink()

    @property
    def assembly_file(self):
        import xcoll as xc
        if self._is_null:
            return None
        if self.is_generic():
            return self._generic_assembly_file
        fedb = xc.fluka.environment.fedb
        file = fedb / "assemblies" / f"{self.fedb_series}_{self.fedb_tag}.lbp"
        return file.resolve()

    @assembly_file.setter
    def assembly_file(self, path):
        import xcoll as xc
        if self._is_null:
            raise ValueError("Cannot set assembly_file for a null assembly!")
        if self.is_defunct():
            raise ValueError(f"Cannot set assembly_file for defunct assembly '{self.name}'!")
        if self.is_generic():
            raise ValueError(f"Cannot set assembly_file for generic assembly '{self.name}'!\n"
                           + f"It is generated automatically.")
        path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        fedb = xc.fluka.environment.fedb
        target = fedb / "assemblies" / f"{self.fedb_series}_{self.fedb_tag}.lbp"
        if path != target:
            path.copy_to(target, method='mount')
        # Store
        dct = self.to_dict()
        # Remove fluka_name for previously auto-generated material,
        # such that the material code can be re-generated
        if 'material' in dct and dct['material'] is not None and 'fluka_name' in dct['material'] \
        and dct['material']['fluka_name'] is not None and dct['material']['fluka_name'].startswith('XCOLL'):
            dct['material'].pop('fluka_name')
        file = fedb / "metadata" / f'{self.fedb_series}_{self.fedb_tag}.lbp.json'
        json.dump(dct, file, indent=2)

    @property
    def body_file(self):
        pass

    @body_file.setter
    def body_file(self, path):
        pass

    @property
    def material_file(self):
        pass

    @material_file.setter
    def material_file(self, path):
        pass

    @property
    def region_file(self):
        pass

    @region_file.setter
    def region_file(self, path):
        pass

    def populate_into_temp_fedb(self, fedb):
        fedb = FsPath(fedb).resolve()
        if not fedb.exists():
            raise ValueError(f"FEDB path {fedb} does not exist!")
        if self.is_generic():
            from xcoll.scattering_routines.fluka.generic_prototype import (
                _assembly_file, _body_file, _crystal_body_file, _region_file,
                _crystal_region_file, _material_file, _crystal_material_file
            )
            if self.is_crystal:
                body_file, tank_file = _crystal_body_file(fedb, self.fedb_tag,
                    self.length, self.bending_radius, self.width, self.height)
                body_region_file, tank_region_file = _crystal_region_file(fedb, self.fedb_tag)
                body_mat_file, tank_mat_file = _crystal_material_file(fedb, self.fedb_tag, self.material)
            else:
                body_file, tank_file = _body_file(fedb, self.fedb_tag, self.length,
                                                  self.width, self.height)
                body_region_file, tank_region_file = _region_file(fedb, self.fedb_tag)
                body_mat_file, tank_mat_file = _material_file(fedb, self.fedb_tag, self.material)
            for pro in self.prototypes:
                if pro.name.endswith('_B'):
                    pro._generic_body_file = body_file
                    pro._generic_region_file = body_region_file
                    pro._generic_material_file = body_mat_file
                elif pro.name.endswith('_T'):
                    pro._generic_body_file = tank_file
                    pro._generic_region_file = tank_region_file
                    pro._generic_material_file = tank_mat_file
                else:
                    raise ValueError(f"Generic assembly prototype '{pro.name}' has invalid name! "
                                   + "Expected to end with '_B' or '_T'.")
            self._generic_assembly_file = _assembly_file(fedb, self.fedb_tag, self.side)
        else:
            link = fedb / 'assemblies' / self.assembly_file.name
            if not link.exists():
                link.symlink_to(self.assembly_file)
        for prot in self.prototypes:
            prot.populate_into_temp_fedb(fedb)

    @property
    def prototypes(self):
        if not hasattr(self, '_prototypes'):
            if not self.assembly_file.exists():
                raise FileNotFoundError(f"Assembly file {self.assembly_file} not found!")
            prototypes = []
            prototype_found = False
            fedb_series = None
            fedb_tag = None
            # Try maximally 5 times to parse the assembly file (in case it's being written by another process)
            for i in range(5):
                success = True
                with self.assembly_file.open('r') as fid:
                    for line in fid:
                        if line.upper().startswith('PROTOTYPE'):
                            prototype_found = True
                            continue
                        if line.upper().startswith('FEDB_SERIES'):
                            if not prototype_found:
                                if i < 4:
                                    # Retry parsing the file
                                    success = False
                                    break
                                raise ValueError("Corrupt assembly file: FEDB_SERIES without PROTOTYPE.")
                            fedb_series = line.split()[1]
                        if line.upper().startswith('FEDB_TAG'):
                            if not prototype_found:
                                if i < 4:
                                    # Retry parsing the file
                                    success = False
                                    break
                                raise ValueError("Corrupt assembly file: FEDB_TAG without PROTOTYPE.")
                            fedb_tag = line.split()[1]
                        if fedb_series is not None and fedb_tag is not None:
                            prototypes.append(FlukaPrototype(fedb_series=fedb_series, fedb_tag=fedb_tag,
                                                            _allow_generic=True))
                            fedb_series = None
                            fedb_tag = None
                            prototype_found = False
                if success:
                    break
                else:
                    time.sleep(0.1)
            if prototype_found:
                raise ValueError("Corrupt assembly file: incomplete prototype definition.")
            self._prototypes = prototypes
        return self._prototypes

    @property
    def files(self):
        if self.assembly_file.exists():
            files = [self.assembly_file]
            for prot in self.prototypes:
                files += prot.files
            return files

    def check_file_valid(self, raise_error=True):
        if self._file_is_valid is None or self._file_is_valid == False:
            if self.assembly_file is None:
                self._file_is_valid = False
                if raise_error:
                    raise ValueError("Assembly has no assembly_file defined!")
                return self._file_is_valid
            self._file_is_valid = True
            prototype_found = False
            fedb_series = None
            fedb_tag = None
            ass_fedb_tag = None
            # Try maximally 5 times to parse the assembly file (in case it's being written by another process)
            for i in range(5):
                success = True
                with self.assembly_file.open('r') as fid:
                    for line in fid:
                        if line.upper().startswith('ASSEMBLY'):
                            if ass_fedb_tag:
                                self._file_is_valid = False
                                if raise_error:
                                    if i < 4:
                                        # Retry parsing the file
                                        success = False
                                        break
                                    raise ValueError("Corrupt assembly file: ASSEMBLY defined more than once.")
                            ass_fedb_tag = line.split()[1]
                            continue
                        if line.upper().startswith('PROTOTYPE'):
                            if ass_fedb_tag:
                                self._file_is_valid = False
                                if raise_error:
                                    if i < 4:
                                        # Retry parsing the file
                                        success = False
                                        break
                                    raise ValueError("Corrupt assembly file: PROTOTYPE after ASSEMBLY.")
                            prototype_found = True
                            continue
                        if line.upper().startswith('FEDB_SERIES'):
                            if fedb_series is not None:
                                self._file_is_valid = False
                                if raise_error:
                                    if i < 4:
                                        # Retry parsing the file
                                        success = False
                                        break
                                    raise ValueError("Corrupt assembly file: FEDB_SERIES defined more than once.")
                            if not prototype_found:
                                self._file_is_valid = False
                                if raise_error:
                                    if i < 4:
                                        # Retry parsing the file
                                        success = False
                                        break
                                    raise ValueError("Corrupt assembly file: FEDB_SERIES without PROTOTYPE.")
                            fedb_series = line.split()[1]
                        if line.upper().startswith('FEDB_TAG'):
                            if fedb_tag is not None:
                                self._file_is_valid = False
                                if raise_error:
                                    if i < 4:
                                        # Retry parsing the file
                                        success = False
                                        break
                                    raise ValueError("Corrupt assembly file: FEDB_TAG defined more than once.")
                            if not prototype_found:
                                self._file_is_valid = False
                                if raise_error:
                                    if i < 4:
                                        # Retry parsing the file
                                        success = False
                                        break
                                    raise ValueError("Corrupt assembly file: FEDB_TAG without PROTOTYPE.")
                            fedb_tag = line.split()[1]
                        if fedb_series is not None and fedb_tag is not None:
                            fedb_series = None
                            fedb_tag = None
                            prototype_found = False
                if not ass_fedb_tag:
                    self._file_is_valid = False
                    if raise_error:
                        if i < 4:
                            # Retry parsing the file
                            success = False
                        else:
                            raise ValueError("Corrupt assembly file: ASSEMBLY not defined.")
                if ass_fedb_tag != self.fedb_tag:
                    self._file_is_valid = False
                    if raise_error:
                        if i < 4:
                            # Retry parsing the file
                            success = False
                        else:
                            raise ValueError(f"Corrupt assembly file: ASSEMBLY {ass_fedb_tag} "
                                    + f"does not match {self.fedb_tag}. Please take note "
                                    + f"that the filename should match exactly, i.e. "
                                    + f"{self.fedb_series}_{self.fedb_tag}.lbp")
                if success:
                    break
                else:
                    time.sleep(0.1)
        return self._file_is_valid


class FlukaPrototypeAccessor:
    """This class is used to access the prototypes in the FEDB."""

    def __init__(self):
        self._type = 'Prototype'

    @property
    def _raw(self):
        return [(pro.fedb_series.lower(), pro.fedb_tag.lower(), pro)
                for pro in FlukaPrototype._registry
                if not isinstance(pro, FlukaAssembly)]

    def __repr__(self):
        return f"<FlukaPrototypeAccessor at {hex(id(self))} (use .show() to see the content)>"

    def __str__(self):
        res = []
        for ss, vv in self.ordered_data.items():
            res.append(f"FEDB Series: {ss}")
            for tt, pro in vv.items():
                if pro.is_broken:
                    continue
                res.append(f"    {tt:<16}: {pro}")
        res_broken = []
        for ss, vv in self.ordered_data.items():
            if any(pro.is_broken for pro in vv.values()):
                res_broken.append(f"FEDB Series: {ss}")
            for tt, pro in vv.items():
                if not pro.is_broken:
                    continue
                res_broken.append(f"    {tt:<16}: {pro}")
        if res_broken:
            res.append('')
            res.append('BROKEN:')
            res.append('')
            res += res_broken
        return "\n".join(res)

    def show(self):
        """Print the content of the FEDB."""
        print(self)

    @property
    def series(self):
        return {ss for ss, _, _ in self._raw}

    @property
    def tags(self):
        return {tt for _, tt, _ in self._raw}

    @property
    def data(self):
        return {f'{ss}_{tt}': pro for ss, tt, pro in self._raw}

    @property
    def ordered_data(self):
        return {
            ss: {tt: vv for xx, tt, vv in self._raw if xx == ss}
            for ss in self.series
        }

    def keys(self):
        f"""A set-like object providing a view on the {self._type.lower()} tags."""
        return self.data.keys()

    def values(self):
        f"""A set-like object providing a view on the {self._type.lower()[-1]}ies."""
        return self.data.values()

    def items(self):
        f"""A set-like object providing a view on the {self._type.lower()[-1]}ies and their tags."""
        return self.data.items()

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return iter(self.values().__iter__())

    def __iter__(self):
        super().__setattr__('_iter_data', iter(self.keys()))
        return self

    def __next__(self):
        try:
            name = next(self._iter_data)
        except StopIteration:
            raise StopIteration
        else:
            return self[name]

    def __contains__(self, val):
        val = val.lower()
        return val in self.data or val in self.series or val in self.tags

    def __getitem__(self, val):
        val = val.lower()
        # If val is the full specification (series_tag), return the prototype
        if val in self.data:
            return self.data[val]
        # If val is a series, return the prototypes in that series
        elif val in self.series or val == 'generic':   # Special case for Andre
            return FlukaSeriesAccessor(self._type, self.ordered_data, val)
        # If val is a tag and that tag is unique, return the prototype
        elif val in self.tags:
            this_data = [pro for ss, tt, pro in self._raw if tt == val]
            if len(this_data) > 1:
                raise KeyError(f"Tag '{val}' is not unique in the FEDB. "
                              + "Please specify the series as well.")
            else:
                return this_data[0]
        else:
            raise KeyError(f"{self._type} '{val}' not found in the FEDB.")


class FlukaAssemblyAccessor(FlukaPrototypeAccessor):
    """This class is used to access the assemblies in the FEDB."""

    def __init__(self):
        self._type = 'Assembly'

    @property
    def _raw(self):
        return [(pro.fedb_series.lower(), pro.fedb_tag.lower(), pro)
                for pro in FlukaPrototype._registry
                if isinstance(pro, FlukaAssembly)]


class FlukaSeriesAccessor:
    """This class is used to access the prototypes or assemblies in a specific series."""

    def __init__(self, type, data, series):
        self._type = type
        self._series = series
        if series == 'generic' and series not in data:
            self._series_data = {}
        else:
            self._series_data = data[series]

    def __repr__(self):
        return f"<FlukaSeriesAccessor at {hex(id(self))} (use .show() to see the content)>"

    def __str__(self):
        res = [f"FEDB Series: {self._series}"]
        for tt, pro in self._series_data.items():
            if pro.is_broken:
                continue
            res.append(f"    {tt:<16}: {pro}")
        res_broken = []
        for tt, pro in self._series_data.items():
            if not pro.is_broken:
                continue
            res_broken.append(f"    {tt:<16}: {pro}")
        if res_broken:
            res.append('')
            res.append('BROKEN:')
            res.append('')
            res += res_broken
        return "\n".join(res)

    def show(self):
        """Print the content of the FEDB series."""
        print(self)

    def keys(self):
        f"""A set-like object providing a view on the {self._type.lower()} tags."""
        return self._series_data.keys()

    def values(self):
        f"""A set-like object providing a view on the {self._type.lower()[-1]}ies."""
        return self._series_data.values()

    def items(self):
        f"""A set-like object providing a view on the {self._type.lower()[-1]}ies and their tags."""
        return self._series_data.items()

    def __len__(self):
        return len(self._series_data)

    def __iter__(self):
        return iter(self.values().__iter__())

    def __contains__(self, val):
        val = val.lower()
        return val in self._series_data

    def __getitem__(self, val):
        val = val.lower()
        # If val is a tag, return the prototype
        if val in self._series_data:
            return self._series_data[val]
        raise KeyError(f"{self._type} with tag '{val}' not found in series "
                     + f"'{self._series}'.")
