# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import io
import xtrack as xt


class BaseWrapper:
    _engine_cls = None
    _environment_cls = None

    def __init__(self):
        self._engine = None
        self._environment = None

    @property
    def engine(self):
        self._lazy_load_engine()
        return self._engine

    @property
    def environment(self):
        self._lazy_load_environment()
        return self._environment

    def _lazy_load_environment(self):
        # Ensure the environment is loaded
        if not self._environment:
            self._environment = _environment_cls()

    def _lazy_load_engine(self):
        # Ensure the engine is loaded
        if not self._engine:
            self._engine = _engine_cls()
            self._engine._environment = self.environment
