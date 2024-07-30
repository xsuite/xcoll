#!/bin/bash

ver=$(gcc -dumpversion)
if [ ${ver%%.*} -lt 9 ]
then
  echo "ERROR: Need GCC 9 or higher!"
  exit 1
fi

# source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos8-gcc11-opt/setup.sh

cd xcoll/scattering_routines/k2/FORTRAN_src
rm *.mod *.o


# Compile libraries
if [[ "$OSTYPE" == "darwin"* ]]
then
  cp crlibm/CMakeLists.txt crlibm/CMakeLists.txt.bak
  cp roundctl/CMakeLists.txt roundctl/CMakeLists.txt.bak
  sed -i '' 's/-mfpmath=sse -msse2//g' crlibm/CMakeLists.txt
  sed -i '' 's/-mfpmath=sse -msse2//g' roundctl/CMakeLists.txt
  CFLAGS="'-fPIC'"
else
  CFLAGS="'-fPIC -O3 -mfpmath=sse -msse2 -mavx -mavx2 -mno-fma4 -mno-fma'"
fi
cd crlibm
make clean
rm -r CMakeCache*.txt CMakeFiles Makefile* cmake_install*.cmake &> /dev/null
cmake .
make CFLAGS=$CFLAGS
mv libcrlibm.a ../
cd ../roundctl
make clean
rm -r CMakeCache*.txt CMakeFiles Makefile* cmake_install*.cmake &> /dev/null
cmake .
make CFLAGS=$CFLAGS
mv libroundctl.a ../
cd ..
if [[ "$OSTYPE" == "darwin"* ]]
then
  mv crlibm/CMakeLists.txt.bak crlibm/CMakeLists.txt
  mv roundctl/CMakeLists.txt.bak roundctl/CMakeLists.txt
fi


# Compile FORTRAN
if [[ "$OSTYPE" == "darwin"* ]]
then
  FFLAGS="-m64 -fpic -funroll-loops -std=f2008 -cpp -DDOUBLE_MATH -DCRLIBM -DROUND_NEAR -O3"
else
  FFLAGS="-m64 -fpic -funroll-loops -std=f2008 -cpp -DDOUBLE_MATH -DCRLIBM -DROUND_NEAR -O3 -mfpmath=sse -msse2 -mavx -mavx2 -mno-fma4 -mno-fma"
fi
gfortran $FFLAGS -c \
 core_tools.f90 \
 constants.f90 \
 strings.f90 \
 mod_alloc.f90 \
 common_modules.f90 \
 string_tools.f90 \
 mod_units.f90 \
 extra.f90 \
 mod_particles.f90 \
 prror.f90 \
 mod_meta.f90 \
 mod_time.f90 \
 bouncy_castle.f90 \
 libcrlibm.a \
 libroundctl.a \
 coll_jawfit.f90 \
 coll_common.f90 \
 coll_db.f90 \
 mod_ranlux.f90 \
 mod_funlux.f90 \
 coll_crystal.f90 \
 coll_k2.f90 \
 coll_dist.f90 \
 collimation.f90


# Link FORTRAN
f2py -m pyk2f -c pyk2.f90 \
 core_tools.o \
 constants.o \
 strings.o \
 mod_alloc.o \
 common_modules.o \
 string_tools.o \
 mod_units.o \
 extra.o \
 mod_particles.o \
 prror.o \
 mod_meta.o \
 mod_time.o \
 bouncy_castle.o \
 libcrlibm.a \
 libroundctl.a \
 coll_jawfit.o \
 coll_common.o \
 coll_db.o \
 mod_ranlux.o \
 mod_funlux.o \
 coll_crystal.o \
 coll_k2.o \
 coll_dist.o \
 collimation.o

mv pyk2f.*.so ../
