# This script needs to be ran only once, to set the FLUKA environment.
import xcoll as xc
from pathlib import Path

# Set the paths to FLUKA executables
xc.fluka.interface.fluka       = '/eos/project/f/flukafiles/fluka-coupling/fluka4-5.0/bin/rfluka'
xc.fluka.interface.flukaserver = '/eos/project/f/flukafiles/fluka-coupling/fluka_coupling/fluka/flukaserver'
xc.fluka.interface.linebuilder = '/eos/project/c/collimation-team/software/linebuilder'
xc.fluka.interface.flair       = '/eos/project/f/flukafiles/fluka-coupling/flair-3.4/flair'

print(xc.fluka.interface)
print()

# Compile the FLUKA interface. This should be done within the environment it will be used to avoid dependency issues.
# E.g. when running on HTCondor with cvmfs, compile with cvmfs sourced.
xc.fluka.interface.compile(flukaio_path='/eos/project-c/collimation-team/software/flukaio', verbose=True)

# Import a FLUKA FEDB
xc.fluka.interface.import_fedb(fedb_path='/eos/project/c/collimation-team/software/fedb_coupling', verbose=True, overwrite=False)
