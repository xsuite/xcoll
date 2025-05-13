# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .reference_masses import fluka_masses
from .reference_names import fluka_names


class FlukaWrapper:
    """Wrapper for all FLUKA functions."""

    def __init__(self):
        self._engine = None
        self._environment = None
        self._masses = fluka_masses
        self._particle_names = fluka_names

    @property
    def engine(self):
        if not self._engine:
            # Lazy load the engine to avoid circular imports
            # and to ensure that the environment is set up first
            from xcoll.scattering_routines.fluka.engine import FlukaEngine
            self._engine = FlukaEngine()
        return self._engine

    @property
    def environment(self):
        if not self._environment:
            # Lazy load the engine to avoid circular imports
            # and to ensure that the environment is set up first
            from xcoll.scattering_routines.fluka.environment import FlukaEnvironment
            self._environment = FlukaEnvironment()
        return self._environment

    @property
    def masses(self):
        return self._masses

    @property
    def particle_names(self):
        return self._particle_names

    @property
    def assemblies(self):
        from xcoll import FlukaPrototype, FlukaAssembly
        filtered = [(pro.fedb_series.lower(), pro.fedb_tag.lower(), pro)
                    for pro in FlukaPrototype._registry
                    if isinstance(pro, FlukaAssembly)]
        # build a dict‐of‐dicts by nesting comprehensions
        series = {k for k, _, _ in filtered}
        return {
            kk: {ii: vv for xx, ii, vv in filtered if xx == kk}
            for kk in series
        }

    @property
    def prototypes(self):
        from xcoll import FlukaPrototype, FlukaAssembly
        filtered = [(pro.fedb_series.lower(), pro.fedb_tag.lower(), pro)
                    for pro in FlukaPrototype._registry
                    if not isinstance(pro, FlukaAssembly)]
        # build a dict‐of‐dicts by nesting comprehensions
        series = {k for k, _, _ in filtered}
        return {
            kk: {ii: vv for xx, ii, vv in filtered if xx == kk}
            for kk in series
        }

    def compile(self, *args, **kwargs):
        """Compile the FLUKA code."""
        self.environment.compile(*args, **kwargs)

    def import_fedb(self, *args, **kwargs):
        """Import a FLUKA FEDB."""
        self.environment.import_fedb(*args, **kwargs)
