# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import io
import xtrack as xt


class BaseWrapper:
    _engine_cls = None
    _environment_cls = None
    _particle_mass_cls = None

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
        if self._particle_mass_cls is None:
            raise NotImplementedError(
                "This wrapper does not implement particle_masses.")
        return self._particle_mass_cls()

    def _lazy_load_environment(self):
        # Ensure the environment is loaded
        if not self._environment:
            self._environment = self._environment_cls()

    def _lazy_load_engine(self):
        # Ensure the engine is loaded
        if not self._engine:
            self._engine = self._engine_cls()
            self._engine._environment = self.environment
            if self._particle_mass_cls is None:
                self._engine._masses = None
            else:
                self._engine._masses = self.particle_masses
