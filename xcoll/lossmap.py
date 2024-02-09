# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #
import numpy as np
import pandas as pd
from pathlib import Path
import json

import xtrack as xt
import xpart as xp
import xobjects as xo

from .beam_elements import _all_collimator_types

# TODO: should we add a check to make sure line_is_reversed is a boolean/correct value?

class LossMap:

    def __init__(self, line, line_is_reversed, part):
        self._line = line
        self._line_is_reversed = line_is_reversed
        self._machine_length = line.get_length()
        self._lossmap = None                            
        self._summary = None
        self._part = part

    # TODO: Check if _all_collimator_types is valid
   
    def lossmap(self, file = None, interpolation = 0.1, weights = None):
        """ 
            Does the loss location refinement and interpolation. 
            Writes the data to file if given one. 
        """
        # loss location refinement
        if interpolation is not None: 
            new_state = self._part.state                   
            new_elem  = self._part.at_element               

            for idx, elem in enumerate(self._part.at_element):
                if (self._part.state[idx] == 0 and elem > 0 
                    and self._line.element_names[elem - 1] in 
                    self._line.get_elements_of_type(_all_collimator_types)[1]):

                    print(f"Found at {self._line.element_names[elem]}, should be {self._line.element_names[elem-1]}")
                    new_elem[idx] = elem - 1
                    what_type = self._line.element_names[elem-1].__class__.__name__
                    if what_type   == 'EverestCollimator':
                        new_state[idx] = -333
                    elif what_type == 'EverestCrystal':   
                        new_state[idx] = -334
                    elif what_type == 'BlackAbsorber':  
                        new_state[idx] = -340
                    else:
                        raise ValueError(f"Unknown collimator type {what_type}")
                    
            self._part.state      = new_state
            self._part.at_element = new_elem

            # interpolation 
            aper_s = list(self._part.s[self._part.state==0])
            if len(aper_s) > 0:
                print("Performing the aperture losses refinement.")
                loss_loc_refinement = xt.LossLocationRefinement(self._line, 
                                                                n_theta = 360, # Angular resolution in the polygonal approximation of the aperture
                                                                r_max = 0.5,   # Maximum transverse aperture [m]
                                                                dr = 50e-6,    # Transverse loss refinement accuracy [m]
                                                                ds = interpolation) # Longitudinal loss refinement accuracy [m]
                loss_loc_refinement.refine_loss_location(self._part)
        aper_s, aper_names, aper_nabs = self._get_aperture_losses(weights)
        coll_summary = self.summary(weights, show_zeros=False).to_dict('list')

        self._lossmap = {
                'collimator': {
                    's':      coll_summary['s'],
                    'name':   coll_summary['collname'],
                    'length': coll_summary['length'],
                    'n':      coll_summary['nabs']
                }
                ,
                'aperture': {
                    's':    aper_s,
                    'name': aper_names,
                    'n':    aper_nabs
                }
                ,
                'machine_length': self._machine_length
                ,
                'interpolation': interpolation
                ,
                'reversed': self._line_is_reversed
            }
        return self._lossmap


    def summary(self, weights = None, show_zeros = True, file = None):
        """
            Returns the summary of the losses in the machine.
            Writes the data to .txt file if given one. 
        """
        collimator_names = self._line.get_elements_of_type(_all_collimator_types)[1] 

        if weights is None:
            weights = np.ones(len(self._part.x))
        else:
            self._part.sort(interleave_lost_particles = True)
        
        coll_mask     = (self._part.state <= -333) & (self._part.state >= -340)
        coll_losses   = np.array([self._line.element_names[i] for i in self._part.at_element[coll_mask]])
        coll_lengths  = [self._line[j].active_length for j in collimator_names] 
        coll_pos      = [(self._line.get_s_position(i) + self._line[i].active_length/2) for i in collimator_names]
        
        if self._line is reversed:
            coll_pos  = [self._machine_length - s for s in coll_pos]
        
        coll_types    = [self._line[i].__class__.__name__ for i in collimator_names]  
        coll_weights  = weights[coll_mask]
        nabs          = [coll_weights[coll_losses == j].sum() for j in collimator_names]

        self._summary = pd.DataFrame({
            'collname': collimator_names, 
            'nabs':   nabs,                   # of particles lost on collimators
            'length':   coll_lengths,
            's':        coll_pos,    
            'type':     coll_types
        })

        if file is not None:
            with open(Path(file), 'w') as fid:
                fid.write(self._summary.__repr__())
        
        if show_zeros:
            return self._summary
        else:
            return self._summary[self._summary.nabs > 0]


    def _get_aperture_losses(self, weights = None):
        """
            Returns the aperture losses.
        """
        if weights is None:
            weights  = np.ones(len(self._part.x))
        else:
            self._part.sort(interleave_lost_particles=True)

        # Get s position per particle (lost on aperture)
        aper_mask    = self._part.state == 0
        aper_s       = list(self._part.s[aper_mask])

        if len(aper_s) == 0:
            return [], [], []
        if self._line_is_reversed:
            aper_s   = [ self._machine_length - s for s in aper_s ]

        # Store names of aperture markers
        aper_names   = [self._line.element_names[i] for i in self._part.at_element[aper_mask]]
        name_dict    = dict(zip(aper_s, aper_names)) # TODO: not floating-point-safe and slow

        # Create output arrays
        aper_pos     = np.unique(aper_s)
        aper_weights = weights[aper_mask]
        aper_nabs    = [aper_weights[aper_s == j].sum() for j in aper_pos] 
        aper_names   = [name_dict[ss] for ss in aper_pos]

        return aper_pos, aper_names, aper_nabs

    def dump(self, file=None):
        """
            Dumps the lossmap to a file.
        """
        if file is not None:
            with open(Path(file), 'w') as fid:
                json.dump(self._lossmap, fid, indent=True, cls=xo.JEncoder)
        else:
            raise ValueError("No file given to dump the lossmap to.")