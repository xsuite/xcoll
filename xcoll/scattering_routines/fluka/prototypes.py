# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from pathlib import Path

from xobjects import HybridClass, String

from .paths import fedb


class FlukaPrototype(HybridClass):
    _xofields = {
        '_fedb_series': String,
        '_fedb_tag':    String,
        '_name':        String,
    }

    # We have a registry for FlukaPrototypes and another for FlukaAssemblies
    _registry = {}

    # A HybridClass needs something to depend on, otherwise the class is added
    # twice in the cdefs during compilation
    _depends_on = [String]

    _store_in_to_dict = ['_info', '_extra_commands']

    def __new__(cls, **kwargs):
        # If the prototype is already active, return the existing instance
        fedb_series = kwargs.get('fedb_series', None)
        fedb_tag = kwargs.get('fedb_tag', None)
        if fedb_series is not None or fedb_tag is not None:
            _registry = {**FlukaPrototype._registry, **FlukaAssembly._registry}
            for prototype in _registry.values():
                if prototype.fedb_series == fedb_series and prototype.fedb_tag == fedb_tag:
                    if not isinstance(prototype, cls):
                        raise ValueError(f"{cls.__name__} '{fedb_series}_{fedb_tag}' already exists "
                                       + f"as a {prototype.__class__.__name__}!")
                    return prototype
        try:
            return HybridClass.__new__(cls, **kwargs)
        except TypeError:
            return HybridClass.__new__(cls)

    def __init__(self, **kwargs):
        to_assign = {}
        if '_xobject' not in kwargs:
            kwargs['_fedb_series'] = kwargs.pop('fedb_series', None)
            kwargs['_fedb_tag'] = kwargs.pop('fedb_tag', None)
            if kwargs['_fedb_series'] is None and kwargs['_fedb_tag'] is None:
                kwargs['_fedb_series'] = ''
                kwargs['_fedb_tag'] = ''
                kwargs['info'] = None
                kwargs['extra_commands'] = None
            elif kwargs['_fedb_series'] is None or kwargs['_fedb_tag'] is None:
                raise ValueError("Both 'fedb_series' and 'fedb_tag' must be provided.")
            kwargs['_fedb_series'] = kwargs['_fedb_series'].ljust(8)  # Pre-allocate 8 chars using whitespace
            kwargs['_fedb_tag'] = kwargs['_fedb_tag'].ljust(8)        # Pre-allocate 8 chars using whitespace
            kwargs['_name'] = kwargs['_fedb_tag']
            to_assign['_info'] = kwargs.pop('info', None)
            to_assign['_extra_commands'] = kwargs.pop('extra_commands', None)
        super().__init__(**kwargs)
        if not hasattr(self, '_id'):
            self._id = None
        if not hasattr(self, '_is_null'):
            if self._fedb_series.strip() == '' and self._fedb_tag.strip() == '':
                self._is_null = True
            else:
                self._is_null = False
        if not hasattr(self, '_type'):
            self._type = self.__class__.__name__[5:].lower()
        if not hasattr(self, '_elements'):
            self._elements = []
        for key, val in to_assign.items():
            setattr(self, key, val)
        # print(self)

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

    # def to_dict(self):
    #     if self._is_null:
    #         return None
    #     return {
    #         '__class__': self.__class__.__name__,
    #         'name': self.name,
    #         'fedb_series': self.fedb_series,
    #         'fedb_tag': self.fedb_tag,
    #         'info': self.info,
    #         'extra_commands': self.extra_commands,
    #     }

    # @classmethod
    # def from_dict(cls, data):
    #     cls = data.pop('__class__', None)
    #     if cls == 'FlukaPrototype':
    #         return FlukaPrototype(**data)
    #     elif cls == 'FlukaAssembly':
    #         return FlukaAssembly(**data)
    #     else:
    #         raise ValueError(f"Invalid data format for {cls}.")

    @property
    def name(self):
        if self._is_null:
            return None
        return self._name.strip()

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
        return self._fedb_series.strip()

    @property
    def fedb_tag(self):
        if self._is_null:
            return None
        return self._fedb_tag.strip()

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
        _registry = {**FlukaPrototype._registry, **FlukaAssembly._registry}
        # Verify that the element is not already assigned to another prototype
        for prototype in _registry.values():
            if prototype is self:
                continue
            if element in prototype.elements:
                raise ValueError(f"Element '{element}' already assigned {prototype.name} "
                               + f"{prototype._type}!")
        # Add the prototype to the registry of active prototypes if not yet present
        if len(self._elements) == 0:
            if not self.exists():
                raise ValueError(f"{self._type.capitalize()} '{self.name}' "
                               + f"does not exist in the FEDB!")
            if self.name in _registry:
                # Rename the prototype if the name is already in use
                i = 0
                while True:
                    new_name = f"{self.name}{i}"
                    if new_name not in _registry:
                        print(f"Warning: {self._type.capitalize()} name {self.name} "
                            + f"already in use. Renaming to '{new_name}'.")
                        self.name = new_name
                        break
                    i += 1
            self._id = self._get_next_id()
            self._registry[self.name] = self
        # Add the element to the list of elements that use this prototype
        if element not in self._elements:
            self._elements.append(element)

    def remove_element(self, element):
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
            # Remove the prototype from the registry of active prototypes
            self._registry.pop(prototype.name)
            # Update the IDs of the remaining prototypes and assemblies
            for this_prototype in FlukaPrototype._registry.values():
                if this_prototype._id > self._id:
                    this_prototype._id -= 1
            for this_prototype in FlukaAssembly._registry.values():
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

    def exists(self):
        if self._is_null:
            return False
        return (fedb / "bodies" / f"{self.fedb_series}_{self.fedb_tag}.bodies").exists()

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


    @classmethod
    def _get_next_id(cls):
        # The IDs should be unique over all prototypes and assemblies
        _registry = {**FlukaPrototype._registry, **FlukaAssembly._registry}
        if len(_registry) == 0:
            return 0
        else:
            return max({prototype._id for prototype in _registry.values()}) + 1

    @classmethod
    def make_prototypes(cls, save=True, path=None):
        prototypes = ["#...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8"]
        # First the prototypes
        for prototype in FlukaPrototype._registry.values():
            if prototype.active:
                assert isinstance(prototype, FlukaPrototype)
                prototypes.append(prototype.generate_code())
        # Then the assemblies
        for prototype in FlukaAssembly._registry.values():
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


class FlukaAssembly(FlukaPrototype):
    # We have a registry for FlukaPrototypes and another for FlukaAssemblies
    _registry = {}

    def exists(self):
        return (fedb / "assemblies" / f"{self.fedb_series}_{self.fedb_tag}.lbp").exists()

# import xcoll as xc
# xc.fluka_assemblies['lhc_tcp'].add_element('TCP.C6L7.B1')
# xc.fluka_assemblies['lhc_tcp'].add_element('TCP.D6L7.B1')
# xc.fluka_assemblies['lhc_tcspm'].add_element('TCSPM.E5R7.B1')
# xc.fluka_assemblies['lhc_tcdqaa'].add_element('TCDQ.A5R7.B1')
# xc.fluka_assemblies['lhc_tcl'].add_element('TCL.4R1.B1')
# xc.fluka_assemblies['lhc_tcl1'].add_element('TCL1.4R1.B1')
# xc.fluka_assemblies['lhc_tcdqab'].add_element('TCDQ.B5R7.B1')
# print(xc.FlukaAssembly.make_prototypes(save=False))

assemblies = {
    # SPS assemblies
    'sps_tcsm': FlukaAssembly(fedb_series='sps', fedb_tag='TCSM',    info="test collimator (hollow jaw)"),
    # LHC assemblies
    'lhc_tcp':    FlukaAssembly(fedb_series='lhc', fedb_tag='TCP',     info="primary with jaw in CFC"),
    'lhc_tcsg':   FlukaAssembly(fedb_series='lhc', fedb_tag='TCSG',    info="secondary with jaw in CFC"),
    'lhc_tcspm':  FlukaAssembly(fedb_series='lhc', fedb_tag='TCSPMC',  info="secondary with jaw in MoGr coated (5um)"),
    'lhc_tcsp':   FlukaAssembly(fedb_series='lhc', fedb_tag='TCSP',    info="secondary with jaw in CFC and in-jaw BPMs (IR6)"),
    'lhc_tcla':   FlukaAssembly(fedb_series='lhc', fedb_tag='TCLA',    info="shower absorber"),
    'lhc_tct':    FlukaAssembly(fedb_series='lhc', fedb_tag='TCT',     info="tertiary"),
    'lhc_tcl':    FlukaAssembly(fedb_series='lhc', fedb_tag='TCL',     info="physics debris absorber"),
    # 'lhc_tcl1':    FlukaAssembly(fedb_series='lhc1', fedb_tag='TCL',     info="physics debris absorber"),
    'lhc_tdi':    FlukaAssembly(fedb_series='lhc', fedb_tag='TDI',     info="injection protection"),
    'lhc_tclia':  FlukaAssembly(fedb_series='lhc', fedb_tag='TCLIA',   info="injection protection"),
    'lhc_tclib':  FlukaAssembly(fedb_series='lhc', fedb_tag='TCLIB',   info="injection protection"),
    'lhc_tcdqaa': FlukaAssembly(fedb_series='lhc', fedb_tag='TCDQnAA', info="dump protection"),
    # 'lhc_tcdqaa': FlukaPrototype(fedb_series='lhc', fedb_tag='TCDQnAA', info="dump protection"),
    'lhc_tcdqab': FlukaAssembly(fedb_series='lhc', fedb_tag='TCDQnAB', info="dump protection"),
    'lhc_tcdqac': FlukaAssembly(fedb_series='lhc', fedb_tag='TCDQnAC', info="dump protection"),
    # HL-LHC assemblies
    'hilumi_tcppm':  FlukaAssembly(fedb_series='hilumi', fedb_tag='TCPPM',    info="primary with jaw in MoGr coated"),
    'hilumi_tcspm':  FlukaAssembly(fedb_series='hilumi', fedb_tag='TCSPM',    info="secondary with jaw in MoGr coated (6um)"),
    'hilumi_tcspmp': FlukaAssembly(fedb_series='hilumi', fedb_tag='TCSPMPRT', info="secondary TCSPM prototype (three stripes)"),
    'hilumi_tcsp':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCSPGRC',  info="secondary with jaw in CFC and in-jaw BPMs (IR6), with Cu coating layer (3um)"),
    'hilumi_tcld':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCLD',     info="shower absorber"),
    'hilumi_tdisp2': FlukaAssembly(fedb_series='hilumi', fedb_tag='TDISP2',   info="injection protection"),
    'hilumi_tdisp8': FlukaAssembly(fedb_series='hilumi', fedb_tag='TDISP8',   info="injection protection"),
    'hilumi_tctx':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCTx',     info="tertiary"),
    'hilumi_tcty':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCTy',     info="tertiary"),
    'hilumi_tclx':   FlukaAssembly(fedb_series='hilumi', fedb_tag='TCLX',     info="physics debris absorber"),
    # FCC assemblies
    'fcc_tcp':  FlukaAssembly(fedb_series='fcc', fedb_tag='TCP',  info="primary"),
    'fcc_tcsg': FlukaAssembly(fedb_series='fcc', fedb_tag='TCSG', info="secondary"),
    'fcc_tcl':  FlukaAssembly(fedb_series='fcc', fedb_tag='TCL',  info="absorber"),
    'fcc_tcdq': FlukaAssembly(fedb_series='fcc', fedb_tag='TCDQ', info="dump protection"),
}
