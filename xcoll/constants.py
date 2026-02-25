# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from .general import _pkg_root
from .headers.particle_states import XcollParticleStates as _XPS
from .interaction_record.interaction_types import XcollInteractions as _XI

_XPS.export_to_module(__name__, include_groups=True, include_src=True,
                      include_meta=True, import_vars=True)
_XI.export_to_module(__name__, include_groups=True, include_src=True,
                     include_meta=True, import_vars=True)

_XPS.export_src(_pkg_root / "generated_src" / "particle_states.h")
_XI.export_src(_pkg_root / "generated_src" / "interaction_types.h")
