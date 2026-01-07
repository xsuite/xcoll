# This script needs to be ran only once, to set the FLUKA environment.
import xcoll as xc
from pathlib import Path

path = Path('/Users/frederik/ExternalLibs')

# Set the paths to FLUKA executables
xc.fluka.environment.fluka	 = path / 'fluka4-5.1' / 'bin' / 'rfluka'
xc.fluka.environment.flukaserver = path / 'fluka4-5.1' / 'bin' / 'flukaserver'
xc.fluka.environment.linebuilder = path / 'linebuilder'
xc.fluka.environment.flair       = path / 'flair-3.4' / 'flair'

print(xc.fluka.environment)
print()

# Compile the FLUKA interface. This should be done within the environment it will be used to avoid dependency issues.
# E.g. when running on HTCondor with cvmfs, compile with cvmfs sourced.
xc.fluka.environment.compile(flukaio_lib=path / 'fluka4-5.1' / 'lib' / 'libFlukaIO.a', verbose=True)

# Import a FLUKA FEDB
xc.fluka.environment.import_fedb(fedb_path=path / 'fedb_coupling', verbose=True, overwrite=False)

# Delete any leftover generic assemblies and prototypes from previous runs
xc.fluka.reset_generic()
