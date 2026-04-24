# Compile for debug
export CFLAGS="-g3 -fsanitize=address -static-libasan "$CFLAGS
export LDFLAGS="-g3 -fsanitize=address -static-libasan "$LDFLAGS
export LD_PRELOAD=/home/fvanderv/miniforge3/envs/xcoll/lib/libasan.so.8
xsuite-prebuild r
gdb pytest test_adt.py::test_blow_up\[B1H-ContextCPU\]

# Only run absorber/Everest/FLUKA/Geant4 tests
pytest -n 12 -m black
pytest -n 12 -m everest
pytest -n 12 -m fluka
pytest -n 12 -m geant4

# Combining marks
pytest -n 12 -m "fluka or geant4"

# Run everything except FLUKA/Geant4
pytest -n 12 -m "not fluka and not geant4"

# Run other xcoll tests
pytest -n 12 -m xcother

# Run xaux tests
pytest -n 12 -m xaux
