#!/bin/bash

ver=$(gcc -dumpversion)
if [ ${ver%%.*} -lt 9 ]
then
  echo "ERROR: Need GCC 9 or higher!"
  exit 1
fi

path=xcoll/scattering_routines/fluka
# source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos8-gcc11-opt/setup.sh

git submodule update --init --recursive

cd $path
rm pyflukaf.*.so

cd flukaio
make libs BUILD64=Y

cd ../FORTRAN_src
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
f2py -m pyflukaf -c pyfluka.f90 \
 core_tools.o \
 constants.o \
 strings.o \
 mod_alloc.o \
 common_modules.o \
 string_tools.o  \
 mod_units.o \
 pdgid.o \
 mod_fluka.o \
 ../flukaio/lib/libFlukaIO64.a

mv pyflukaf.* ../
cd ../../../..
echo
if [ -f ${path}/pyflukaf.*.so ]
then
    echo "Created pyFLUKA shared library in "$( ls ${path}/pyflukaf.*.so )
else
    echo "Failed pyFLUKA compilation! No shared library found in "${path}" !"
fi
