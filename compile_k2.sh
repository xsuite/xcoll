#!/bin/bash

ver=$(gcc -dumpversion)
if [ ${ver%%.*} -lt 9 ]
then
  echo "ERROR: Need GCC 9 or higher!"
  exit 1
fi

if [[ "$OSTYPE" == "darwin"* ]]
then
  echo "ERROR: Cannot compile K2 on MacOS."
  exit 1
fi


# source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos8-gcc11-opt/setup.sh

cd xcoll/scattering_routines/k2/FORTRAN_src
rm *.mod *.o

# Compile libraries
CFLAGS='-fPIC -O3 -mfpmath=sse -msse2 -mavx -mavx2 -mno-fma4 -mno-fma'
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

# Compile FORTRAN
files="core_tools.f90 constants.f90 strings.f90 mod_alloc.f90 common_modules.f90 string_tools.f90 mod_units.f90 extra.f90 mod_particles.f90 prror.f90 mod_meta.f90 mod_time.f90 bouncy_castle.f90 libcrlibm.a libroundctl.a coll_jawfit.f90 coll_common.f90 coll_db.f90 mod_ranlux.f90 mod_funlux.f90 coll_crystal.f90 coll_k2.f90 coll_dist.f90 collimation.f90"
gfortran -m64 -fpic -funroll-loops -std=f2008 -cpp -DDOUBLE_MATH -DCRLIBM -DROUND_NEAR -O3 -mfpmath=sse -msse2 -mavx -mavx2 -mno-fma4 -mno-fma -c $files

# Link FORTRAN
f2py -m pyk2f -c pyk2.f90 ${files/.f90/.o}

mv pyk2f.*.so ../
