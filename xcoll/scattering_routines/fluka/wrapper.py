# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import io
import xtrack as xt
try:
    from xaux import FsPath  # TODO: once xaux is in Xsuite keep only this
except (ImportError, ModuleNotFoundError):
    from ...xaux import FsPath

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

    def view(self, elements=None, *, input_file=None):
        if elements:
            if input_file:
                raise ValueError("Cannot view elements and input file at the same time!")
            part =xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=1e9)
            input_file = self.engine.generate_input_file(elements=elements, particle_ref=part)
        if not hasattr(input_file, '__iter__') or isinstance(input_file, (str, io.IOBase)):
            input_file = [input_file]
        self.environment.run_flair(input_file[0])
        if elements:
            # Input file was generated from elements, so we can remove it
            for f in input_file:
                if FsPath(f).exists():
                    f.unlink()

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
