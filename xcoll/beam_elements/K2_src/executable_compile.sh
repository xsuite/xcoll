#!/bin/bash
rm *.mod *.o

gfortran \
core_tools.f90 constants.f90 strings.f90 mod_alloc.f90 common_modules.f90 string_tools.f90 mod_units.f90 bouncy_castle.f90  coll_jawfit.f90 coll_common.f90 coll_db.f90 mod_ranlux.f90 mod_funlux.f90 coll_crystal.f90 coll_k2.f90 \
files.f90 main.f90 libcrlibm.a libroundctl.a -g -O0 -o main