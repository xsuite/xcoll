#!/bin/bash

ver=$(gcc -dumpversion)
if [ ${ver%%.*} -lt 9 ]
then
  echo "ERROR: Need GCC 9 or higher!"
  exit 1
fi

# source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos8-gcc11-opt/setup.sh

cd xcoll/scattering_routines/fluka/FORTRAN_src

rm *.mod *.o

gfortran -fpic -c \
 core_tools.f90 \
 constants.f90 \
 strings.f90 \
 mod_alloc.f90 \
 common_modules.f90 \
 string_tools.f90  \
 mod_units.f90 \
 pdgid.f90 \
 mod_fluka.f90

# link fortran
#f2py -m pyFLUKAf -c pyFLUKA.f90 \
# core_tools.o \
# constants.o \
# strings.o \
# mod_alloc.o \
# common_modules.o  \
# string_tools.o  \
# mod_units.o  \
# bouncy_castle.o  \
# libcrlibm.a  \
# libroundctl.a  \
# coll_jawfit.o  \
# coll_common.o  \
# mod_ranlux.o  \
# mod_funlux.o  \
# coll_crystal.o  \
# coll_k2.o \

cd ../../../..
