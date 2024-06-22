# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2024.                 #
# ######################################### #

from pathlib import Path
import numpy as np


def create_dat_file(line, names, file="k2_colldb.dat"):
    from ...beam_elements import _K2Collimator, _K2Crystal
    file = Path(file).expanduser().resolve()
    if file.exists():
        file.unlink()
    with file.open('w') as fp:
        onesided = []
        crystal  = []
        for name in names:
            el = line[name]
            if el.__class__ == _K2Crystal:
                orbit = el.co[0]
            elif el.__class__ == _K2Collimator:
                orbit = el.co[0][0]
            else:
                print(f"Warning: ignoring {name} in K2 input file generation as it "
                      f"is not an _K2Collimator or _K2Crystal.")
                continue
            if not el.optics_ready():
                raise ValueError(f"Collimator {name} has no optics assigned! "
                                 f"Do this before the K2 input file generation.")
            if hasattr(el.angle, '__iter__'):
                raise ValueError(f"Collimator {name} has jaw-dependent angles. "
                                 f"Not supported.")
            gap = el.gap
            if not isinstance(gap, str) and hasattr(gap, '__iter__'):
                gap = gap[0]  # Only left-sided collimators supported in SixTrack
            if gap == None:
                print(f"Warning: ignoring {name} in K2 input file generation as it "
                      f"has no assigned opening.")
                continue
            mat = el.material
            length = el.length
            angle = el.angle
            # SixTrack offset is wrt. closed orbit
            if el.side == 'both':
                offset = np.round((el.jaw_L + el.jaw_R)/2 - orbit, 9) + 0 # adding 0 to avoid -0 output
            elif el.side == 'left':
                offset = 0.
                onesided.append(name)
            else:
                raise ValueError(f"Right-sided collimator {name} not supported.")
            if el.__class__ == _K2Crystal:
                # Have to do this here and not at the start, to avoid having appended the crystal list
                # but then continuing (e.g. because the gap is None)
                crystal.append(name)
            fp.write(f"{name:50} {gap:>15} {mat:10} {length:<8} {angle:>5}   {offset}\n")

        if len(onesided) > 0:
            fp.write("SETTINGS \n")
            for name in onesided:
                el = line[name]
                fp.write(f"ONESIDED  {name:50}     {el._side}\n")

        if len(crystal) > 0:
            if len(onesided) == 0:
                fp.write("SETTINGS \n")
            for name in crystal:
                el = line[name]
                fp.write(f"CRYSTAL   {name:18}  {el.bending_radius:5}  {el.width}  "
                         f"{el.height}  0.0  {el.tilt}  {el.miscut}  {el._orient}\n")

    print(f"Created SixTrack CollDB {file}")
    return file
