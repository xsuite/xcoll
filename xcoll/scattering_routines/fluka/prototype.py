# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import json
import numpy as np
from pathlib import Path

try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

from .environment import format_fluka_float
from ...beam_elements.base import BaseCollimator


class FlukaPrototype:
    # This is a registry to keep track of prototypes. If one is already defined,
    # we do not create a new instance but return the existing one
    _registry = []

    # This is a registry to keep track of prototypes that are in use (and for which
    # we will need to generate protoype code). We have a separate registry for
    # FlukaPrototypes and another for FlukaAssemblies.
    _assigned_registry = {}

    def __new__(cls, **kwargs):
        # If the prototype is already defined, return the existing instance
        fedb_series = kwargs.get('fedb_series', None)
        fedb_tag = kwargs.get('fedb_tag', None)
        if fedb_series is not None or fedb_tag is not None:
            for prototype in FlukaPrototype._registry:
                if prototype.fedb_series.upper() == fedb_series.upper() \
                and prototype.fedb_tag.upper() == fedb_tag.upper():
                    if prototype.__class__ != cls:
                        raise ValueError(f"{cls.__name__} '{fedb_tag}' is already defined "
                                       + f"as a {prototype.__class__.__name__}!")
                    return prototype
        # Register the new prototype
        self = object.__new__(cls)
        if fedb_series is not None and fedb_tag is not None:
            FlukaPrototype._registry.append(self)
        return self

    def __init__(self, fedb_series=None, fedb_tag=None, *, angle=0, side=None, width=None,
                 height=None, length=None, material=None, info=None, extra_commands=None,
                 is_crystal=False, bending_radius=None, **kwargs):
        if fedb_series is None and fedb_tag is None:
            self._is_null = True
            info = None
            extra_commands = None
        elif fedb_series is None or fedb_tag is None:
            raise ValueError("Both 'fedb_series' and 'fedb_tag' must be provided.")
        else:
            self._is_null = False
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
        self._material = material  # This will overwrite the material in the FlukaCollimator
        self._is_crystal = is_crystal
        self._bending_radius = bending_radius
        self._info = info
        self._extra_commands = extra_commands
        self._id = None
        self._type = self.__class__.__name__[5:].lower()
        self._elements = []

    def __repr__(self):
        if self._is_null:
            return ''
        if self.assigned:
            n_active = len(self.active_elements)
            n_total = len(self.elements)
            n_inactive = n_total - n_active
            elements = []
            if n_active > 0:
                elements.append(f"{n_active} active")
            if n_inactive > 0:
                elements.append(f"{n_inactive} inactive")
            elements = " and ".join(elements)
            elements += " element"
            if n_total > 1:
                elements += 's'
        else:
            elements = "unassigned"
        info = f" ({self.info})" if self.info else ''
        return f"{self.__class__.__name__} '{self.name}' ({elements}): " \
             + f"tag {self.fedb_tag} in {self.fedb_series} series{info}"

    def __str__(self):
        if self._is_null:
            return ''
        return self.__repr__()

    def delete(self):
        if self._is_null:
            return
        if len(self._elements) > 0:
            raise ValueError(f"Cannot delete {self._type} '{self.name}' "
                           + f"while it has {len(self._elements)} elements assigned!")
        # Remove the prototype from the registry of all prototypes
        if self in FlukaPrototype._registry:
            FlukaPrototype._registry.remove(self)
        # Remove the prototype from the registry of all assemblies
        if self in FlukaAssembly._registry:
            FlukaAssembly._registry.remove(self)
        # Remove all files associated with the prototype
        for file in self.files:
            file.unlink(missing_ok=True)

    def to_dict(self):
        if self._is_null:
            return {'__class__': self.__class__.__name__}
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
            'material': self.material,
            'is_crystal': self.is_crystal,
            'bending_radius': self.bending_radius,
            'info': self.info,
            'extra_commands': self.extra_commands,
        }

    @classmethod
    def from_dict(cls, data):
        cls = data.pop('__class__', None)
        if cls == 'FlukaPrototype':
            return FlukaPrototype(**data)
        elif cls == 'FlukaAssembly':
            return FlukaAssembly(**data)
        elif cls == 'FlukaGenericAssembly':
            from xcoll import FlukaGenericAssembly
            return FlukaGenericAssembly(**data)
        elif cls == 'FlukaGenericCrystalAssembly':
            from xcoll import FlukaGenericCrystalAssembly
            return FlukaGenericCrystalAssembly(**data)
        else:
            raise ValueError(f"Invalid data format for {cls}.")

    @classmethod
    def from_json(cls, path):
        if isinstance(path, str):
            path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        with path.open('r') as fid:
            data = json.load(fid)
        return cls.from_dict(data)

    @property
    def name(self):
        if self._is_null:
            return None
        return self._name

    @name.setter
    def name(self, val):
        if self._is_null:
            return
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
        if self._is_null:
            return None
        import xcoll as xc
        file = xc.fluka.environment.fedb / "bodies" \
                    / f"{self.fedb_series}_{self.fedb_tag}.bodies"
        return file.resolve()

    @body_file.setter
    def body_file(self, path):
        if self._is_null:
            raise ValueError("Cannot set body_file for a null prototype!")
        import xcoll as xc
        path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        path.copy_to(xc.fluka.environment.fedb / "bodies" / f"{self.fedb_series}_{self.fedb_tag}.bodies")
        with open(xc.fluka.environment.fedb / "metadata" / f'{self.fedb_series}_{self.fedb_tag}.bodies.json', 'w') as fid:
            json.dump(self.to_dict(), fid, indent=4)

    @property
    def material_file(self):
        if self._is_null:
            return None
        import xcoll as xc
        file = xc.fluka.environment.fedb / "materials" \
                    / f"{self.fedb_series}_{self.fedb_tag}.assignmat"
        return file.resolve()

    @material_file.setter
    def material_file(self, path):
        if self._is_null:
            raise ValueError("Cannot set material_file for a null prototype!")
        import xcoll as xc
        path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        path.copy_to(xc.fluka.environment.fedb / "materials" / f"{self.fedb_series}_{self.fedb_tag}.assignmat")

    @property
    def region_file(self):
        if self._is_null:
            return None
        import xcoll as xc
        file = xc.fluka.environment.fedb / "regions" \
                    / f"{self.fedb_series}_{self.fedb_tag}.regions"
        return file.resolve()

    @region_file.setter
    def region_file(self, path):
        if self._is_null:
            raise ValueError("Cannot set region_file for a null prototype!")
        import xcoll as xc
        path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        path.copy_to(xc.fluka.environment.fedb / "regions" / f"{self.fedb_series}_{self.fedb_tag}.regions")

    @property
    def files(self):
        return [self.body_file, self.material_file, self.region_file]

    def exists(self):
        if self._is_null:
            return False
        return np.all([ff.exists() for ff in self.files])

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
    def assigned(self):
        if self._is_null:
            return False
        return self._id is not None

    @property
    def active(self):
        if self._is_null:
            return False
        return len(self.active_elements) > 0

    @property
    def fluka_position(self):
        # Maximum positions for the parking region:
        # x in [-3000.0, 3000.0], y in [-4000.0, -2000.0], z in [0.0, 1.E5]
        if self.assigned:
            return 0., 0., 0., (self._id%5-2)*500.0 , -3000., (self._id//5+1)*1000.0
        else:
            return None

    @property
    def elements(self):
        return self._elements.copy()

    @property
    def active_elements(self):
        active_elements = []
        for ee in self.elements:
            if (not hasattr(ee, 'active') or ee.active) \
            and (not hasattr(ee, 'jaw') or ee.jaw):
                active_elements.append(ee)
        return active_elements

    def add_element(self, element, force=True):
        from .engine import FlukaEngine
        if element is None:
            if not force:
                raise ValueError("Cannot add a null element to a prototype!")
            return None
        elif not hasattr(element, 'name'):
            if not force:
                raise ValueError(f"Element {element} has no `name` variable! "
                               + f"Cannot assign to a prototype!")
            return None
        elif not isinstance(element, FlukaEngine()._element_classes):
            if not force:
                raise ValueError(f"Element {element.name} is not a FLUKA element! "
                               + f"Cannot assign to a prototype!")
            return None
        if self._is_null:
            if not force:
                raise ValueError(f"Cannot add element {element.name} to a null prototype!")
            return None
        _registry = {**FlukaPrototype._assigned_registry, **FlukaAssembly._assigned_registry}
        # Verify that the element is not already assigned to another prototype
        for prototype in _registry.values():
            if prototype is self:
                continue
            if element in prototype.elements:
                raise ValueError(f"Element '{element}' already assigned {prototype.name} "
                               + f"{prototype._type}!")
        if not self.exists():
            raise ValueError(f"{self._type.capitalize()} '{self.name}' "
                           + f"does not exist in the FEDB!")
        # Add the prototype to the registry of assigned prototypes if not yet present
        if len(self._elements) == 0:
            # Rename the prototype if the name is already in use
            existing_names = [key.upper() for key in _registry.keys()]
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
            self._id = self._get_next_id()
            self._assigned_registry[self.name] = self
        # Add the element to the list of elements that use this prototype
        if element not in self._elements:
            self._elements.append(element)

    def remove_element(self, element, force=True):
        if self._is_null:
            if force:
                raise ValueError("Cannot remove element from a null prototype!")
            return None
        if element is None:
            if force:
                raise ValueError("Cannot remove a null element from a prototype!")
            return None
        # Remove the element from the list of elements that use this prototype
        if element in self._elements:
            self._elements.remove(element)
        if len(self._elements) == 0:
            if self.name in self._assigned_registry:
                # Remove the prototype from the registry of active prototypes
                self._assigned_registry.pop(self.name)
                # Update the IDs of the remaining prototypes and assemblies
                for this_prototype in FlukaPrototype._assigned_registry.values():
                    if this_prototype._id > self._id:
                        this_prototype._id -= 1
                for this_prototype in FlukaAssembly._assigned_registry.values():
                    if this_prototype._id > self._id:
                        this_prototype._id -= 1
                self._id = None

    def generate_code(self):
        if self.active and self.active_elements:
            _type = 'ASSEMBLY' if isinstance(self, FlukaAssembly) else 'PROTOTYPE'
            prot  = f"{_type:9}     {self.name}\n"
            prot += f"FEDB_SERIES   {self.fedb_series}\n"
            prot += f"FEDB_TAG      {self.fedb_tag}\n"
            prot += f"ROT-DEFI  "
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
            for element in self.active_elements:
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
        else:
            return '#'


    def view(self, show=True, keep_files=False):
        import xcoll as xc
        xc.fluka.environment.test_assembly(self.fedb_series, self.fedb_tag, show=show,
                                           keep_files=keep_files)


    @classmethod
    def _get_next_id(cls):
        # The IDs should be unique over all prototypes and assemblies
        _assigned_registry = {**FlukaPrototype._assigned_registry, **FlukaAssembly._assigned_registry}
        if len(_assigned_registry) == 0:
            return 0
        else:
            return max({prototype._id for prototype in _assigned_registry.values()}) + 1

    @classmethod
    def make_prototypes(cls, save=True, path=None):
        prototypes = ["#...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8"]
        # First the prototypes
        for prototype in FlukaPrototype._assigned_registry.values():
            if prototype.active:
                assert isinstance(prototype, FlukaPrototype)
                prototypes.append(prototype.generate_code())
        # Then the assemblies
        for prototype in FlukaAssembly._assigned_registry.values():
            if prototype.active:
                assert isinstance(prototype, FlukaAssembly)
                prototypes.append(prototype.generate_code())
        prototypes = "\n".join(prototypes)
        if save:
            if path is None:
                path = Path.cwd()
            with (path / "prototypes.lbp").open("w") as fp:
                fp.write(prototypes)
        return prototypes

    @classmethod
    def inspect_prototypes_file(cls, prototypes_file):
        if prototypes_file is None:
            prototypes_file = Path.cwd() / "prototypes.lbp"
        all_prototypes = {}
        with Path(prototypes_file).open("r") as fp:
            name = None
            for line in fp.readlines():
                line = line.strip()
                if line.startswith("#") or len(line) == 0:
                    continue
                if line.startswith("PROTOTYPE"):
                    name = line.split()[1]
                    all_prototypes[name] = {'type': 'PROTOTYPE'}
                elif line.startswith("ASSEMBLY"):
                    name = line.split()[1]
                    all_prototypes[name] = {'type': 'ASSEMBLY'}
                elif line.startswith("FEDB_SERIES"):
                    if name is None:
                        raise ValueError(f"File {prototypes_file} malformed ("
                                    + "FEDB_SERIES without PROTOTYPE or ASSEMBLY)!")
                    all_prototypes[name]['fedb_series'] = line.split()[1]
                elif line.startswith("FEDB_TAG"):
                    if name is None:
                        raise ValueError(f"File {prototypes_file} malformed ("
                                    + "FEDB_TAG without PROTOTYPE or ASSEMBLY)!")
                    all_prototypes[name]['fedb_tag'] = line.split()[1]
                elif line.startswith("MAP_ENTRIES"):
                    if name is None:
                        raise ValueError(f"File {prototypes_file} malformed ("
                                    + "MAP_ENTRIES without PROTOTYPE or ASSEMBLY)!")
                    if 'elements' not in all_prototypes[name]:
                        all_prototypes[name]['elements'] = []
                    all_prototypes[name]['elements'] += line.split()[1:]
        for prototype in [*FlukaPrototype._assigned_registry.values(),
                          *FlukaAssembly._assigned_registry.values()]:
            if prototype.name not in all_prototypes:
                if prototype.active_elements:
                    raise ValueError(f"Prototype {prototype.name} not found in prototypes file!")
            else:
                _type = 'ASSEMBLY' if isinstance(prototype, FlukaAssembly) else 'PROTOTYPE'
                if _type != all_prototypes[prototype.name]['type']:
                    raise ValueError(f"Wrong type for {prototype.name} (expected "
                                  + f"{prototype._type.upper()}, got "
                                  + f"{all_prototypes[prototype.name]['type']})!")
                if 'fedb_series' not in all_prototypes[prototype.name]:
                    raise ValueError(f"FEDB_SERIES for {prototype.name} not found in prototypes file!")
                if prototype.fedb_series != all_prototypes[prototype.name]['fedb_series']:
                    raise ValueError(f"Wrong fedb_series for {prototype.name} (expected "
                                   + f"{prototype.fedb_series}, got "
                                   + f"{all_prototypes[prototype.name]['fedb_series']})!")
                if 'fedb_tag' not in all_prototypes[prototype.name]:
                    raise ValueError(f"FEDB_TAG for {prototype.name} not found in prototypes file!")
                if prototype.fedb_tag != all_prototypes[prototype.name]['fedb_tag']:
                    raise ValueError(f"Wrong fedb_tag for {prototype.name} (expected "
                                   + f"{prototype.fedb_tag}, got "
                                   + f"{all_prototypes[prototype.name]['fedb_tag']})!")
                if 'elements' not in all_prototypes[prototype.name]:
                    raise ValueError(f"MAP_ENTRIES for {prototype.name} not found in prototypes file!")
                for element in prototype.active_elements:
                    if element.name.upper() not in all_prototypes[prototype.name]['elements']:
                        raise ValueError(f"Element {element.name} not found in prototypes file!")


class FlukaAssembly(FlukaPrototype):
    # We have a registry for FlukaPrototypes and another for FlukaAssemblies
    _assigned_registry = {}

    @property
    def assembly_file(self):
        import xcoll as xc
        file = xc.fluka.environment.fedb / "assemblies" \
                    / f"{self.fedb_series}_{self.fedb_tag}.lbp"
        return file.resolve()

    @assembly_file.setter
    def assembly_file(self, path):
        if self._is_null:
            raise ValueError("Cannot set assembly_file for a null prototype!")
        import xcoll as xc
        path = FsPath(path)
        if not path.exists():
            raise FileNotFoundError(f"File {path} does not exist!")
        path.copy_to(xc.fluka.environment.fedb / "assemblies" / f"{self.fedb_series}_{self.fedb_tag}.lbp")
        with open(xc.fluka.environment.fedb / "metadata" / f'{self.fedb_series}_{self.fedb_tag}.lbp.json', 'w') as fid:
            json.dump(self.to_dict(), fid, indent=4)

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

    @property
    def prototypes(self):
        if not self.assembly_file.exists():
            raise FileNotFoundError(f"Assembly file {self.assembly_file} not found!")
        prototypes = []
        prototype_found = False
        fedb_series = None
        fedb_tag = None
        with self.assembly_file.open('r') as fid:
            for line in fid:
                if line.upper().startswith('PROTOTYPE'):
                    prototype_found = True
                    continue
                if line.upper().startswith('FEDB_SERIES'):
                    if not prototype_found:
                        raise ValueError("Corrupt assembly file: FEDB_SERIES without PROTOTYPE.")
                    fedb_series = line.split()[1]
                if line.upper().startswith('FEDB_TAG'):
                    if not prototype_found:
                        raise ValueError("Corrupt assembly file: FEDB_TAG without PROTOTYPE.")
                    fedb_tag = line.split()[1]
                if fedb_series is not None and fedb_tag is not None:
                    prototypes.append(FlukaPrototype(fedb_series=fedb_series, fedb_tag=fedb_tag))
                    fedb_series = None
                    fedb_tag = None
                    prototype_found = False
        if prototype_found:
            raise ValueError("Corrupt assembly file: incomplete prototype definition.")
        return prototypes

    @property
    def files(self):
        files = [self.assembly_file]
        for prot in self.prototypes:
            files += prot.files
        return files
