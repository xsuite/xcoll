# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xtrack.particles.pdg as pdg

class BaseWrapper:
    _engine_cls = None
    _environment_cls = None
    _particle_masses_meta = None
    _particle_names_meta = None

    def __init__(self):
        self._engine = None
        self._environment = None

    def __str__(self):
        return f"{self.__class__.__name__}"

    def __repr__(self):
        return f"<{str(self)} at {hex(id(self))}>"

    @property
    def engine(self):
        self._lazy_load_engine()
        return self._engine

    @property
    def environment(self):
        self._lazy_load_environment()
        return self._environment

    @property
    def particle_masses(self):
        if self._particle_masses_meta is None:
            raise NotImplementedError(
                "This wrapper does not implement particle_masses.")
        return BaseWrapperAccessor(self._particle_masses_meta)

    @property
    def particle_names(self):
        if self._particle_names_meta is None:
            raise NotImplementedError(
                "This wrapper does not implement particle_names.")
        return BaseWrapperAccessor(self._particle_names_meta)

    def _lazy_load_environment(self):
        # Ensure the environment is loaded
        if not self._environment:
            self._environment = self._environment_cls()

    def _lazy_load_engine(self):
        # Ensure the engine is loaded
        if not self._engine:
            self._engine = self._engine_cls()
            self._engine._environment = self.environment
            if self._particle_masses_meta is None:
                self._engine._masses = None
            else:
                self._engine._masses = self.particle_masses


class BaseWrapperAccessor:
    """Accessor to get reference masses."""

    def __init__(self, meta_dict):
        self._meta_dict = meta_dict

    def __getitem__(self, pdgid: int | str) -> float | None:
        if isinstance(pdgid, str):
            pdgid = pdg.get_pdg_id_from_name(pdgid)
        vals = [vv['value'] for _, vv in self._meta_dict.items()
                if vv['pdg_id'] == pdgid]
        if len(vals) == 1:
            return vals[0]
        vals = [vv['value'] for _, vv in self._meta_dict.items()
                if vv['pdg_id'] == -pdgid]
        if len(vals) == 1:
            return vals[0]

    def __contains__(self, pdgid: int | str) -> bool:
        return self[pdgid] is not None