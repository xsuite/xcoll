#!/bin/bash

ver=$(gcc -dumpversion)
if [ ${ver%%.*} -lt 9 ]
then
  echo "ERROR: Need GCC 9 or higher!"
  exit 1
fi

# source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos8-gcc11-opt/setup.sh

# Only for testing:
cd xcoll/scattering_routines/pyeverest/FORTRAN_src
rm *.mod *.o

# compile fortran
gfortran -fpic -c \
 mod_ranlux.f90  \
 mod_funlux.f90  \
 coll_k2.f90 \

# link fortran
f2py -m pyk2f -c pyk2.f90 \
 mod_ranlux.o  \
 mod_funlux.o  \
 coll_k2.o \

mv pyk2f.*.so ../

cd ../../../..