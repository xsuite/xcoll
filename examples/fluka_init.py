# This script needs to be ran only once, to set the FLUKA environment.
import xcoll as xc
from pathlib import Path

path = Path('/Users/frederik/ExternalLibs')

# Set the paths to FLUKA executables
xc.fluka.interface.fluka	 = path / 'fluka4-5.1' / 'bin' / 'rfluka'
xc.fluka.interface.flukaserver = path / 'fluka4-5.1' / 'bin' / 'flukaserver'
xc.fluka.interface.linebuilder = path / 'linebuilder'
xc.fluka.interface.flair       = path / 'flair-3.4' / 'flair'

print(xc.fluka.interface)
print()

# Compile the FLUKA interface. This should be done within the environment it will be used to avoid dependency issues.
# E.g. when running on HTCondor with cvmfs, compile with cvmfs sourced.
xc.fluka.interface.compile(flukaio_lib=path / 'fluka4-5.1' / 'lib' / 'libFlukaIO.a', verbose=True)

# Import a FLUKA FEDB
xc.fluka.interface.import_fedb(fedb_path=path / 'fedb_coupling', verbose=True, overwrite=True)
