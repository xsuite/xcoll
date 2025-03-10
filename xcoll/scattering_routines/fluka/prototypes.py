# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import numpy as np
from pathlib import Path

from .environment import FlukaEnvironment
from ...beam_elements.base import BaseCollimator


class FlukaPrototype:
    # This is a registry to keep track of prototypes. If one is already defined,
    # we do not create a new instance but return the existing one
    _registry = []

    # This is a registry to keep track of prototypes that are in use (and for which
    # we will need to generate protoype code). We have a separate registry for
    # FlukaPrototypes and another for FlukaAssemblies.
    _active_registry = {}

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
        FlukaPrototype._registry.append(self)
        return self

    def __init__(self, fedb_series=None, fedb_tag=None, angle=0, side=None, material=None,
                 info=None, extra_commands=None):
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
        self._angle = angle
        self._side = side
        self._material = material
        self._info = info
        self._extra_commands = extra_commands
        self._id = None
        self._type = self.__class__.__name__[5:].lower()
        self._elements = []

    def __repr__(self):
        if self._is_null:
            return None
        if self.active:
            if len(self.elements) == 1:
                elements = f"{len(self.elements)} element"
            else:
                elements = f"{len(self.elements)} elements"
        else:
            elements = "unused"
        return f"{self.__class__.__name__} '{self.name}' ({elements}): " \
             + f"tag {self.fedb_tag} in {self.fedb_series} series ({self.info})"

    def __str__(self):
        if self._is_null:
            return ''
        return self.__repr__()

    def to_dict(self):
        if self._is_null:
            return {'__class__': self.__class__.__name__}
        return {
            '__class__': self.__class__.__name__,
            'name': self.name,
            'fedb_series': self.fedb_series,
            'fedb_tag': self.fedb_tag,
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
        else:
            raise ValueError(f"Invalid data format for {cls}.")

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
        file = FlukaEnvironment().fedb_base / "bodies" \
                    / f"{self.fedb_series}_{self.fedb_tag}.bodies"
        return file.resolve()

    @property
    def material_file(self):
        file = FlukaEnvironment().fedb_base / "materials" \
                    / f"{self.fedb_series}_{self.fedb_tag}.assignmat"
        return file.resolve()

    @property
    def region_file(self):
        file = FlukaEnvironment().fedb_base / "regions" \
                    / f"{self.fedb_series}_{self.fedb_tag}.regions"
        return file.resolve()

    @property
    def files(self):
        return [self.body_file, self.material_file, self.region_file]

    def exists(self):
        if self._is_null:
            return False
        return np.all([ff.exists() for ff in self.files])

    @property
    def angle(self):
        if self._is_null:
            return None
        return self._angle

    @property
    def side(self):
        if self._is_null:
            return None
        return BaseCollimator.side.fget(self)

    @side.setter
    def side(self, val):
        BaseCollimator.side.fset(self)

    @property
    def material(self):
        if self._is_null:
            return None
        return self._material

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
    def active(self):
        if self._is_null:
            return False
        return self._id is not None

    @property
    def fluka_position(self):
        # Maximum positions for the parking region:
        # x in [-3000.0, 3000.0], y in [-4000.0, -2000.0], z in [0.0, 1.E5]
        if self.active:
            return 0., 0., 0., (self._id%5-2)*500.0 , -3000., (self._id//5+1)*1000.0
        else:
            return None

    @property
    def elements(self):
        return self._elements

    def add_element(self, element, force=True):
        if self._is_null:
            if force:
                raise ValueError("Cannot add element to a null prototype!")
            return None
        if element is None:
            if force:
                raise ValueError("Cannot add a null element to a prototype!")
            return None
        _registry = {**FlukaPrototype._active_registry, **FlukaAssembly._active_registry}
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
        # Add the prototype to the registry of active prototypes if not yet present
        if len(self._elements) == 0:
            if not self.exists():
                raise ValueError(f"{self._type.capitalize()} '{self.name}' "
                               + f"does not exist in the FEDB!")
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
            self._active_registry[self.name] = self
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
            if self.name in self._active_registry:
                # Remove the prototype from the registry of active prototypes
                self._active_registry.pop(self.name)
                # Update the IDs of the remaining prototypes and assemblies
                for this_prototype in FlukaPrototype._active_registry.values():
                    if this_prototype._id > self._id:
                        this_prototype._id -= 1
                for this_prototype in FlukaAssembly._active_registry.values():
                    if this_prototype._id > self._id:
                        this_prototype._id -= 1
                self._id = None

    def generate_code(self):
        if self.active:
            prot  = f"{self._type.upper():9}     {self.name}\n"
            prot += f"FEDB_SERIES   {self.fedb_series}\n"
            prot += f"FEDB_TAG      {self.fedb_tag}\n"
            prot += f"ROT-DEFI  "
            for value in self.fluka_position:
                prot += f"{value:>10.1f}"
            prot += " proto\n"
            if self._extra_commands:
                if hasattr(self._extra_commands, "__iter__") \
                and not isinstance(self._extra_commands, str):
                    prot += "\n".join(self._extra_commands) + "\n"
                else:
                    prot += self._extra_commands + "\n"
            maps = []
            for element in self.elements:
                maps.append(f" {element.upper()}")
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

    def in_file(self, file=None):
        if self._is_null:
            return False
        if file is None:
            file = Path.cwd() / "prototypes.lbp"
        file = Path(file)
        if self.name is None or self.fedb_series is None or self.fedb_tag is None:
            return False
        with file.open("r") as fp:
            found_name = False
            found_series = False
            found_tag = False
            for line in fp.readlines():
                if self.name in line and self._type.upper() in line:
                    found_name = True
                if self.fedb_series in line and "FEDB_SERIES" in line:
                    found_series = True
                if self.fedb_tag in line and "FEDB_TAG" in line:
                    found_tag = True
        return found_name and found_series and found_tag


    def view(self, show=True, keep_files=False):
        FlukaEnvironment.test_assembly(self.fedb_series, self.fedb_tag, show=show,
                                       keep_files=keep_files)


    @classmethod
    def _get_next_id(cls):
        # The IDs should be unique over all prototypes and assemblies
        _active_registry = {**FlukaPrototype._active_registry, **FlukaAssembly._active_registry}
        if len(_active_registry) == 0:
            return 0
        else:
            return max({prototype._id for prototype in _active_registry.values()}) + 1

    @classmethod
    def make_prototypes(cls, save=True, path=None):
        prototypes = ["#...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8"]
        # First the prototypes
        for prototype in FlukaPrototype._active_registry.values():
            if prototype.active:
                assert isinstance(prototype, FlukaPrototype)
                prototypes.append(prototype.generate_code())
        # Then the assemblies
        for prototype in FlukaAssembly._active_registry.values():
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
    def reset(cls):
        # First the prototypes
        for prototype in list(FlukaPrototype._active_registry.values()):
            if prototype.active:
                for element in prototype.elements:
                    prototype.remove_element(element)
        assert FlukaPrototype._active_registry == {}
        # Then the assemblies
        for prototype in list(FlukaAssembly._active_registry.values()):
            if prototype.active:
                for element in prototype.elements:
                    prototype.remove_element(element)
        assert FlukaAssembly._active_registry == {}


class FlukaAssembly(FlukaPrototype):
    # We have a registry for FlukaPrototypes and another for FlukaAssemblies
    _active_registry = {}

    @property
    def assembly_file(self):
        file = FlukaEnvironment().fedb_base / "assemblies" \
                    / f"{self.fedb_series}_{self.fedb_tag}.lbp"
        return file.resolve()

    @property
    def body_file(self):
        pass

    @property
    def material_file(self):
        pass

    @property
    def region_file(self):
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


assemblies = {
    # SPS assemblies
    'sps_tcsm':      FlukaAssembly(fedb_series='sps',    fedb_tag='TCSM',     info="test collimator (hollow jaw)"),
    # LHC assemblies
    'lhc_tcp':       FlukaAssembly(fedb_series='lhc',    fedb_tag='TCP',      info="primary with jaw in CFC"),
    'lhc_tcsg':      FlukaAssembly(fedb_series='lhc',    fedb_tag='TCSG',     info="secondary with jaw in CFC"),
    'lhc_tcsp':      FlukaAssembly(fedb_series='lhc',    fedb_tag='TCSP',     info="secondary with jaw in CFC and in-jaw BPMs (IR6)"),
    'lhc_tcla':      FlukaAssembly(fedb_series='lhc',    fedb_tag='TCLA',     info="shower absorber"),
    'lhc_tct':       FlukaAssembly(fedb_series='lhc',    fedb_tag='TCT',      info="tertiary"),
    'lhc_tcl':       FlukaAssembly(fedb_series='lhc',    fedb_tag='TCL',      info="physics debris absorber"),
    'lhc_tdi':       FlukaAssembly(fedb_series='lhc',    fedb_tag='TDI',      angle=np.deg2rad(90.), info="injection protection"),
    'lhc_tclia':     FlukaAssembly(fedb_series='lhc',    fedb_tag='TCLIA',    info="injection protection"),
    'lhc_tclib':     FlukaAssembly(fedb_series='lhc',    fedb_tag='TCLIB',    info="injection protection"),
    'lhc_tcdqaa':    FlukaAssembly(fedb_series='lhc',    fedb_tag='TCDQnAA',  side='left',  info="dump protection"),
    'lhc_tcdqab':    FlukaAssembly(fedb_series='lhc',    fedb_tag='TCDQnAB',  side='left',  info="dump protection"),
    'lhc_tcdqac':    FlukaAssembly(fedb_series='lhc',    fedb_tag='TCDQnAC',  side='left',  info="dump protection"),
    # HL-LHC assemblies
    'hilumi_tcppm':  FlukaAssembly(fedb_series='hilumi', fedb_tag='TCPPM',    info="primary with jaw in MoGr coated"),
    'hilumi_tcspm':  FlukaAssembly(fedb_series='hilumi', fedb_tag='TCSPM',    info="secondary with jaw in MoGr coated (6um)"),
    'hilumi_tcsg':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCSPGRC',  info="secondary with jaw in CFC with Cu coating layer (3um)"),
    'hilumi_tcld':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCLD',     info="shower absorber"),
    'hilumi_tctx':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCTx',     info="tertiary"),
    'hilumi_tcty':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCTy',     info="tertiary"),
    'hilumi_tclx':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCLX',     info="physics debris absorber"),
    # 'hilumi_tcpc':   FlukaPrototype(fedb_series='hilumi',fedb_tag='TCPCHB1',  info="crystal"),
    # FCC assemblies
    'fcc_tcp':       FlukaAssembly(fedb_series='fcc',    fedb_tag='TCP',      info="primary"),
    'fcc_tcsg':      FlukaAssembly(fedb_series='fcc',    fedb_tag='TCSG',     info="secondary"),
    'fcc_tcdq':      FlukaAssembly(fedb_series='fcc',    fedb_tag='TCDQ',     side='right',  info="dump protection"),
    # Only for testing
    'test_donadon':  FlukaAssembly(fedb_series='test',   fedb_tag='DONADON',  info="mesh to test child particle ids"),
    # Old assemblies - not to be used
    'lhc_tcpm':      FlukaAssembly(fedb_series='lhc',    fedb_tag='TCPM',     info="primary with jaw in MoGr coated ??"),
    'lhc_tcspm':     FlukaAssembly(fedb_series='lhc',    fedb_tag='TCSPM',    info="secondary with jaw in MoGr coated"),
    'lhc_tcspmc':    FlukaAssembly(fedb_series='lhc',    fedb_tag='TCSPMC',   info="secondary with jaw in MoGr coated (5um)"),
    'lhc_tcspmp':    FlukaAssembly(fedb_series='lhc',    fedb_tag='TCSPMP',   info="secondary with jaw in MoGr coated (prototype)"),
    'hilumi_tctpx':  FlukaAssembly(fedb_series='hilumi', fedb_tag='TCTPX',    info="alternative tertiary"),
    'hilumi_tctpxv': FlukaAssembly(fedb_series='hilumi', fedb_tag='TCTPXV',   angle=np.deg2rad(90), info="alternative tertiary"),
}

# The following assemblies give wrong results with the jaw test:
assemblies_wrong_jaw = {
    'lhc_tcdqaa_':   FlukaAssembly(fedb_series='lhc',    fedb_tag='TCDQAA',   info="dump protection"),
    'lhc_tcdqab_':   FlukaAssembly(fedb_series='lhc',    fedb_tag='TCDQAB',   info="dump protection"),
    'lhc_tcdqac_':   FlukaAssembly(fedb_series='lhc',    fedb_tag='TCDQAC',   info="dump protection"),
}

# The following assemblies have errors in the prototype code:
assemblies_invalid_prototype = {
    'hilumi_tcspmp': FlukaAssembly(fedb_series='hilumi', fedb_tag='TCSPMPRT', info="secondary TCSPM prototype (three stripes)"),
    'hilumi_tdisp2': FlukaAssembly(fedb_series='hilumi', fedb_tag='TDISP2',   info="injection protection"),
    'hilumi_tdisp8': FlukaAssembly(fedb_series='hilumi', fedb_tag='TDISP8',   info="injection protection"),
    'fcc_tcl':       FlukaAssembly(fedb_series='fcc',    fedb_tag='TCL',      info="absorber"),
}
