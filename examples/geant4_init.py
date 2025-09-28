# This script needs to be ran only once, to set the Geant4 environment.
import xcoll as xc

print(xc.geant4.environment)
xc.geant4.environment.compile()
