# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2023.                 #
# ######################################### #





class LossMap:
    
    def __init__(self, line, part, line_is_reversed):
        self._line = line
        self._line_is_reversed = line_is_reversed
        self._machine_length = line.get_length()
        self._lossmap = None
        self._summary = None



    def _get_aperture_losses(self, part, weights=None):
        if weights is None:
            weights = np.ones(len(part.x))
        else:
            part.sort(interleave_lost_particles=True)

        # Get s position per particle (lost on aperture)
        aper_mask = part.state==0
        aper_s = list(part.s[aper_mask])
        if len(aper_s) == 0:
            return [], [], []
        if self._line_is_reversed:
            aper_s = [ self.machine_length - s for s in aper_s ]

        # Store names of aperture markers
        aper_names   = [self.line.element_names[i] for i in part.at_element[aper_mask]]
        name_dict    = dict(zip(aper_s, aper_names)) # TODO: not floating-point-safe and slow

        # Create output arrays
        aper_pos     = np.unique(aper_s)
        aper_weights = weights[aper_mask]
        aper_nabs    = [aper_weights[aper_s==ss].sum() for ss in aper_pos] # TODO: this might be slow
        aper_names   = [name_dict[ss] for ss in aper_pos]
        return aper_pos, aper_names, aper_nabs