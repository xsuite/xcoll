# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .reference_masses import fluka_masses
from .reference_names import fluka_names
from .prototype import FlukaPrototypeAccessor, FlukaAssemblyAccessor
from .engine import FlukaEngine
from .environment import FlukaEnvironment


class FlukaWrapper:
    """Wrapper for all FLUKA functions."""

    def __init__(self):
        self._engine = None
        self._environment = None
        self._masses = fluka_masses
        self._particle_names = fluka_names

    @property
    def engine(self):
        self._lazy_load_engine()
        return self._engine

    @property
    def environment(self):
        self._lazy_load_environment()
        return self._environment

    @property
    def masses(self):
        return self._masses

    @property
    def particle_names(self):
        return self._particle_names

    @property
    def assemblies(self):
        self._lazy_load_environment()
        return FlukaAssemblyAccessor()

    @property
    def prototypes(self):
        self._lazy_load_environment()
        return FlukaPrototypeAccessor()

    def _lazy_load_environment(self):
        """Ensure the environment is loaded."""
        if not self._environment:
            self._environment = FlukaEnvironment()

    def _lazy_load_engine(self):
        """Ensure the engine is loaded."""
        if not self._engine:
            self._engine = FlukaEngine()
            self._engine._environment = self.environment
        return self._engine
