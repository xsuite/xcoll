#!/bin/bash

ver=$(gcc -dumpversion)
if [ ${ver%%.*} -lt 9 ]
then
  echo "ERROR: Need GCC 9 or higher!"
  exit 1
fi

# source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos8-gcc11-opt/setup.sh

cd scattering_routines/k2/FORTRAN_src
rm *.mod *.o


# compile libraries
cd crlibm
make clean
rm -r CMakeCache*.txt CMakeFiles Makefile* cmake_install*.cmake &> /dev/null
cmake .
make CFLAGS=-fPIC
mv libcrlibm.a ../
cd ../roundctl
make clean
rm -r CMakeCache*.txt CMakeFiles Makefile* cmake_install*.cmake &> /dev/null
cmake .
make CFLAGS=-fPIC
mv libroundctl.a ../
cd ..

# compile fortran
gfortran -fpic -c \
 core_tools.f90 \
 constants.f90 \
 strings.f90 \
 mod_alloc.f90 \
 common_modules.f90  \
 string_tools.f90  \
 mod_units.f90  \
 bouncy_castle.f90  \
 libcrlibm.a \
 libroundctl.a \
 coll_jawfit.f90  \
 coll_common.f90  \
 mod_ranlux.f90  \
 mod_funlux.f90  \
 coll_crystal.f90  \
 coll_k2.f90 \

# link fortran
f2py -m pyk2f -c pyk2.f90 \
 core_tools.o \
 constants.o \
 strings.o \
 mod_alloc.o \
 common_modules.o  \
 string_tools.o  \
 mod_units.o  \
 bouncy_castle.o  \
 libcrlibm.a  \
 libroundctl.a  \
 coll_jawfit.o  \
 coll_common.o  \
 mod_ranlux.o  \
 mod_funlux.o  \
 coll_crystal.o  \
 coll_k2.o \

mv pyk2f.*.so ../
