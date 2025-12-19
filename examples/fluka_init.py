# This script needs to be ran only once, to set the FLUKA environment.
import xcoll as xc
from pathlib import Path

fluka_path = Path('/home/fvanderv/pythondev/fluka4-5.1')
flair_path = Path('/home/fvanderv/pythondev/flair-3.4')

# Set the paths to FLUKA executables
xc.fluka.environment.fluka	 = fluka_path / 'bin' / 'rfluka'
xc.fluka.environment.flukaserver = fluka_path / 'bin' / 'flukaserver'
xc.fluka.environment.linebuilder = '/eos/project/c/collimation-team/software/fluka_coupling_tmp_patch_xsuite'
xc.fluka.environment.flair       = flair_path

print(xc.fluka.environment)
print()

# Compile the FLUKA interface. This should be done within the environment it will be used to avoid dependency issues.
# E.g. when running on HTCondor with cvmfs, compile with cvmfs sourced.
xc.fluka.environment.compile(flukaio_lib=fluka_path / 'interface' / 'flukaio' / 'lib' / 'libFlukaIO64.a', verbose=True)

# Import a FLUKA FEDB
xc.fluka.environment.import_fedb(fedb_path='/eos/project/c/collimation-team/software/fedb_coupling', verbose=True, overwrite=False)
