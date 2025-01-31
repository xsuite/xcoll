# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from pathlib import Path

from .paths import fedb

class FlukaAssembly:
    _registry = {}

    @classmethod
    def get_active_assemblies(cls):
        _ids = cls._registry.copy()
        _ids.pop(-1)
        return _ids

    @classmethod
    def _get_next_id(cls):
        if len(cls._registry) == 0:
            return 0
        else:
            return max({assembly._id for assembly in cls._registry.values()}) + 1

    @classmethod
    def _add_assembly(cls, assembly):
        if assembly.name in cls._registry:
            raise ValueError(f"Assembly '{assembly.name}' already active: {assembly}.")
        if not assembly.exists():
            raise ValueError(f"Assembly '{assembly.name}' does not exist in the FEDB!")
        assembly._id = cls._get_next_id()
        cls._registry[assembly.name] = assembly

    @classmethod
    def _remove_assembly(cls, assembly):
        cls._registry.pop(assembly.name)
        for this_assembly in cls._registry.values():
            if this_assembly._id > assembly._id:
                this_assembly._id -= 1
        assembly._id = None

    @classmethod
    def make_prototypes(cls, save=True, path=Path.cwd()):
        prototypes = []
        for assembly in cls._registry.values():
            if assembly.active:
                prototypes.append(assembly.prototype)
        prototypes = "\n".join(prototypes)
        if save:
            with (path / "prototypes.lbp").open("w") as fp:
                fp.write(prototypes)
        return prototypes

    def __init__(self, name, fedb_series, fedb_tag, info=None):
        self._name = name
        self._fedb_series = fedb_series
        self._fedb_tag = fedb_tag
        self._info = info
        self._entries = []
        self._id = None

    def __repr__(self):
        if self.active:
            if len(self.entries) == 1:
                entries = f"{len(self.entries)} entry"
            else:
                entries = f"{len(self.entries)} entries"
        else:
            entries = "unused"
        return f"FlukaAssembly '{self.name}' ({entries}): " \
             + f"{self.fedb_tag} in {self.fedb_series} ({self.info})"

    def __str__(self):
        return self.__repr__()

    @property
    def name(self):
        return self._name

    @property
    def fedb_series(self):
        return self._fedb_series

    @property
    def fedb_tag(self):
        return self._fedb_tag

    @property
    def info(self):
        return self._info

    @property
    def entries(self):
        return self._entries

    def add_entry(self, entry):
        if len(self._entries) == 0:
            FlukaAssembly._add_assembly(self)
        if entry not in self._entries:
            self._entries.append(entry)

    def remove_entry(self, entry):
        if entry in self._entries:
            self._entries.remove(entry)
        if len(self._entries) == 0:
            FlukaAssembly._remove_assembly(self)

    @property
    def active(self):
        return self._id is not None

    @property
    def fluka_position(self):
        if self.active:
            return 0., 0., 0., (self._id%5-2)*500.0 , -3000., (self._id//5+1)*1000.0
        else:
            return None

    @property
    def prototype(self):
        if self.active:
            prot  = f"ASSEMBLY      {self.name.upper()}\n"
            prot += f"FEDB_SERIES   {self.fedb_series}\n"
            prot += f"FEDB_TAG      {self.fedb_tag}\n"
            prot += f"ROT-DEFI  "
            for value in self.fluka_position:
                prot += f"{value:>10.1f}"
            prot += " proto\n"
            maps = []
            for entry in self.entries:
                maps.append(f" {entry.upper()}")
            if len(maps) > 0:
                prot_map = [f"MAP_ENTRIES  "]
                for entry in maps:
                    if len(prot_map[-1] + entry) > 75:
                        prot_map.append(f"MAP_ENTRIES  ")
                    prot_map[-1] += entry
                prot += "\n".join(prot_map) + "\n"
            prot += "#\n"
            return prot
        else:
            return None

    def exists(self):
        return (fedb / "assemblies" / f"{self.fedb_series}_{self.fedb_tag}.lbp").exists()


assemblies = {
    # LHC assemblies
    'tcp':    FlukaAssembly(name='tcp',    fedb_series='lhc', fedb_tag='TCP',     info="primary with jaw in CFC"),
    'tcpm':   FlukaAssembly(name='tcpm',   fedb_series='lhc', fedb_tag='TCPM',    info="primary with jaw in MoGr coated"),
    'tcsg':   FlukaAssembly(name='tcsg',   fedb_series='lhc', fedb_tag='TCSG',    info="secondary with jaw in CFC"),
    # TODO: is the prototype superflous?
    'tcspmp': FlukaAssembly(name='tcspmp', fedb_series='lhc', fedb_tag='TCSPMC',  info="secondary TCSPM prototype (three stripes)"),
    'tcspm':  FlukaAssembly(name='tcspm',  fedb_series='lhc', fedb_tag='TCSPMC',  info="secondary with jaw in MoGr coated"),
    'tcsp':   FlukaAssembly(name='tcsp',   fedb_series='lhc', fedb_tag='TCSP',    info="secondary of TCSPs in IR6"),
    'tcla':   FlukaAssembly(name='tcla',   fedb_series='lhc', fedb_tag='TCLA',    info="shower absorber"),
    # TODO: can we use one assembly for all TCTs?
    'tctpv':  FlukaAssembly(name='tctpv',  fedb_series='lhc', fedb_tag='TCT',     info="tertiary"),
    'tctph':  FlukaAssembly(name='tctph',  fedb_series='lhc', fedb_tag='TCT',     info="tertiary"),
    'tcl':    FlukaAssembly(name='tcl',    fedb_series='lhc', fedb_tag='TCL',     info="physics debris absorber"),
    'tclia':  FlukaAssembly(name='tclia',  fedb_series='lhc', fedb_tag='TCLIA',   info="injection protection"),
    'tclib':  FlukaAssembly(name='tclib',  fedb_series='lhc', fedb_tag='TCLIB',   info="injection protection"),
    'tdi':    FlukaAssembly(name='tdi',    fedb_series='lhc', fedb_tag='TDI',     info="injection protection"),
    'tcdqaa': FlukaAssembly(name='tcdqaa', fedb_series='lhc', fedb_tag='TCDQnAA', info="extraction protection"),
    'tcdqab': FlukaAssembly(name='tcdqab', fedb_series='lhc', fedb_tag='TCDQnAB', info="extraction protection"),
    'tcdqac': FlukaAssembly(name='tcdqac', fedb_series='lhc', fedb_tag='TCDQnAC', info="extraction protection"),
    # HL-LHC assemblies
    'tcld':   FlukaAssembly(name='tcld',   fedb_series='hilumi', fedb_tag='TCLD', info="shower absorber"),
}
