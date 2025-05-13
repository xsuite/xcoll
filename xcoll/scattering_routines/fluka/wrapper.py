# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .engine import FlukaEngine
from .environment import FlukaEnvironment
from .prototype import FlukaPrototype
from .reference_masses import fluka_masses


class FlukaWrapper:
    """Wrapper for all FLUKA functions."""

    def __init__(self):
        self._engine = FlukaEngine()
        self._environment = FlukaEnvironment()
        self._masses = fluka_masses

    @property
    def engine(self):
        return self._engine

    @property
    def environment(self):
        return self._environment

    @property
    def masses(self):
        return self._masses

    @property
    def assemblies(self):
        return {f'{pro.fedb_series}_{pro.fedb_tag}': pro
                for pro in FlukaPrototype._registry}

    def compile(self, *args, **kwargs):
        """Compile the FLUKA code."""
        self.environment.compile(*args, **kwargs)


fluka = FlukaWrapper()
