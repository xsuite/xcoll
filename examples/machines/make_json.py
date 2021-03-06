#!/usr/bin/env python
# coding: utf-8

# # Prepare HL-LHC Optics

import sys
import json
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xpart as xp

from cpymad.madx import Madx

# Run the mask file to produce the input for xtrack

def run_mask(mask, outfile, sequence):
    mad = Madx()

    mad.call(mask)

    madsequence = getattr(mad.sequence, sequence)

    line = xt.Line.from_madx_sequence(madsequence, apply_madx_errors=True, install_apertures=True)
    line.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, gamma0=madsequence.beam.gamma)

    # Save to json
    with open(outfile, 'w') as fid:
        json.dump(line.to_dict(), fid, cls=xo.JEncoder)

    # Load from json to check
    with open(outfile, 'r') as fid:
        loaded_dct = json.load(fid)
    line_2 = xt.Line.from_dict(loaded_dct)


def main():
    mask = 'lhc_run3_b1.madx'
    outfile = 'lhc_run3_b1.json'
    run_mask(mask, outfile, sequence='lhcb1')
    mask = 'lhc_run3_b2.madx'
    outfile = 'lhc_run3_b2.json'
    run_mask(mask, outfile, sequence='lhcb2')

if __name__ == "__main__":
    main()
