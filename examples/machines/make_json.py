#!/usr/bin/env python

import json
import xobjects as xo
import xtrack as xt
import xpart as xp
from cpymad.madx import Madx

# Run the mask file to produce the input for xtrack

def run_mask(mask, outfile, sequence):
    mad = Madx()
    mad.call(mask)

    line = xt.Line.from_madx_sequence(mad.sequence[sequence], apply_madx_errors=True, install_apertures=True)
    line.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, gamma0=mad.sequence[sequence].beam.gamma)
    
    # Fix cavity frequencies due to bug in mad loader
    cavities, cav_names = line.get_elements_of_type(xt.Cavity)
    for cc, nn in zip(cavities, cav_names):
        if cc.frequency ==0.:
            ii_mad = mad.sequence[sequence].element_names().index(nn)
            cc_mad = mad.sequence[sequence].elements[ii_mad]
            f0_mad = mad.sequence[sequence].beam.freq0 * 1e6 # mad has it in MHz
            harmon = 35640 # cc_mad.parent.harmon
            cc.frequency = f0_mad*harmon

    # Save to json
    with open(outfile, 'w') as fid:
        json.dump(line.to_dict(), fid, cls=xo.JEncoder, indent=True)

    # Load from json to check that there are no loading errors
    with open(outfile, 'r') as fid:
        loaded_dct = json.load(fid)
    line_2 = xt.Line.from_dict(loaded_dct)


def main():
    run_mask(mask='lhc_run3_b1.madx', outfile='lhc_run3_b1_from_mask.json', sequence='lhcb1')
    run_mask(mask='lhc_run3_b2.madx', outfile='lhc_run3_b2_from_mask.json', sequence='lhcb2')

if __name__ == "__main__":
    main()

