# This script needs to be ran only once, to set the FLUKA environment.
import xcoll as xc

# Set the paths to FLUKA executables
xc.fluka.environment.fluka       = '/eos/project/f/flukafiles/fluka-coupling/fluka4-5.0/bin/rfluka'
xc.fluka.environment.flukaserver = '/eos/project/f/flukafiles/fluka-coupling/fluka_coupling/fluka/flukaserver'
xc.fluka.environment.linebuilder = '/eos/project/c/collimation-team/software/fluka_coupling_tmp_patch_xsuite'
xc.fluka.environment.flair       = '/eos/project/f/flukafiles/fluka-coupling/flair-3.4/flair'

# Compile the FLUKA interface
xc.fluka.environment.compile(flukaio_path='/eos/project-c/collimation-team/software/flukaio', verbose=True)

# Import a FLUKA FEDB
xc.fluka.environment.import_fedb(fedb_path='/eos/project/c/collimation-team/software/fedb_coupling', verbose=True, overwrite=False)
