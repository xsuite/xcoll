# copyright ############################### #
# This file is part of the Xcoll package.   #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

import xobjects as xo

from ..general import _pkg_root
from .base import InvalidXcoll


class ChannellingDev(InvalidXcoll):
    _xofields = {
        'length': xo.Float64,
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','channelling.h')
    ]


class BentChannellingDev(InvalidXcoll):
    _xofields = {
        'length': xo.Float64,
        'method' : xo.Int64,  # 2, 3, 4
        'variant': xo.Int64,  # 1 or 2
        'order'  : xo.Int64,  # 2,4,6,8,10,12
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channelling.h')
    ]
               

# Model 2
# --- Order 2 ---
class BentChannellingDevM2V1o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o02.h')
    ]
    
class BentChannellingDevM2V2o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o02.h')
    ] 
    
# --- Order 4 ---
class BentChannellingDevM2V1o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o04.h')
    ]
    
class BentChannellingDevM2V2o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o04.h')
    ] 
    
    
# --- Order 6 ---
class BentChannellingDevM2V1o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o06.h')
    ]
    
class BentChannellingDevM2V2o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o06.h')
    ] 
    
# --- Order 8 ---
class BentChannellingDevM2V1o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o08.h')
    ]
    
class BentChannellingDevM2V2o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o08.h')
    ]

# --- Order 10 ---
class BentChannellingDevM2V1o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o10.h')
    ]
    
class BentChannellingDevM2V2o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o10.h')
    ]

# --- Order 12 ---
class BentChannellingDevM2V1o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V1o12.h')
    ]
    
class BentChannellingDevM2V2o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM2V2o12.h')
    ]


# Model 3
# --- Order 2 ---
class BentChannellingDevM3V1o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o02.h')
    ]
    
class BentChannellingDevM3V2o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o02.h')
    ] 
    
# --- Order 4 ---
class BentChannellingDevM3V1o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o04.h')
    ]
    
class BentChannellingDevM3V2o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o04.h')
    ] 
    
    
# --- Order 6 ---
class BentChannellingDevM3V1o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o06.h')
    ]
    
class BentChannellingDevM3V2o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o06.h')
    ] 
    
# --- Order 8 ---
class BentChannellingDevM3V1o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o08.h')
    ]
    
class BentChannellingDevM3V2o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o08.h')
    ]

# --- Order 10 ---
class BentChannellingDevM3V1o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o10.h')
    ]
    
class BentChannellingDevM3V2o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o10.h')
    ]

# --- Order 12 ---
class BentChannellingDevM3V1o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V1o12.h')
    ]
    
class BentChannellingDevM3V2o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM3V2o12.h')
    ]

# Model 4
# --- Order 2 ---
class BentChannellingDevM4V1o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o02.h')
    ]
    
class BentChannellingDevM4V2o02(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o02.h')
    ] 
    
# --- Order 4 ---
class BentChannellingDevM4V1o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o04.h')
    ]
    
class BentChannellingDevM4V2o04(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o04.h')
    ] 
    
    
# --- Order 6 ---
class BentChannellingDevM4V1o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o06.h')
    ]
    
class BentChannellingDevM4V2o06(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o06.h')
    ] 
    
# --- Order 8 ---
class BentChannellingDevM4V1o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o08.h')
    ]
    
class BentChannellingDevM4V2o08(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o08.h')
    ]

# --- Order 10 ---
class BentChannellingDevM4V1o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o10.h')
    ]
    
class BentChannellingDevM4V2o10(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o10.h')
    ]

# --- Order 12 ---
class BentChannellingDevM4V1o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V1o12.h')
    ]
    
class BentChannellingDevM4V2o12(InvalidXcoll):
    _xofields = {
        'length': xo.Float64
    }

    isthick = True
    behaves_like_drift = False
    allow_track = True
    needs_rng = False
    skip_in_loss_location_refinement = True
    allow_loss_refinement = False

    _depends_on = [InvalidXcoll]

    _extra_c_sources = [
        _pkg_root.joinpath('beam_elements','elements_src','elliptic_functions.h'),
        _pkg_root.joinpath('beam_elements','elements_src','bent_channellingM4V2o12.h')
    ]

                                       
