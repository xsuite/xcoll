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

from .reference_masses import fluka_masses_meta
from .reference_names import fluka_names_meta
from .prototype import FlukaPrototypeAccessor, FlukaAssemblyAccessor
from .engine import FlukaEngine
from .environment import FlukaInterface
from ..wrapper import BaseWrapper


class FlukaWrapper(BaseWrapper):
    """Wrapper for all FLUKA functions."""

    _engine_cls = FlukaEngine
    _interface_cls = FlukaInterface
    _particle_masses_meta = fluka_masses_meta
    _particle_names_meta = fluka_names_meta

    @property
    def assemblies(self):
        self._lazy_load_interface()
        return FlukaAssemblyAccessor()

    @property
    def prototypes(self):
        self._lazy_load_interface()
        return FlukaPrototypeAccessor()

    def view(self, elements=None, *, input_file=None):
        if elements:
            if input_file:
                raise ValueError("Cannot view elements and input file at the same time!")
            part = xt.Particles.reference_from_pdg_id(pdg_id='proton', p0c=1e9)
            input_file = self.engine.generate_input_file(elements=elements, particle_ref=part)
        elif input_file is None:
            raise ValueError("Either elements or input_file must be provided!")
        if not hasattr(input_file, '__iter__') or isinstance(input_file, (str, io.IOBase)):
            input_file = [input_file]
        self.interface.run_flair(input_file[0])
        if elements:
            # Input file was generated from elements, so we can remove it
            for f in input_file:
                if FsPath(f).exists():
                    FsPath(f).unlink()
