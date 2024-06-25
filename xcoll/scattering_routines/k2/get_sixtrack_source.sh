#!/bin/bash

if [ $# -ne 1 ]
then
  echo "This script needs exactly 1 argument: the path to the original SixTrack source."
  exit 1
fi

path=${1%/}/

cd FORTRAN_src

for f in core_tools constants strings mod_alloc common_modules string_tools mod_units mod_particles mod_meta mod_time bouncy_castle coll_jawfit coll_common coll_db mod_ranlux mod_funlux coll_crystal coll_k2 coll_dist collimation
do
    file=${path}${f}.f90
    if [ ! -f $file ]
    then
        echo "Error! Source file $file not found."
        exit 1
    fi
    cp $file .
done

mv crlibm/CMakeLists.txt crlibm_CMakeLists.txt
for file in ${path}crlibm/*
do
    cp $file crlibm/
done
mv crlibm_CMakeLists.txt crlibm/CMakeLists.txt

mv roundctl/CMakeLists.txt roundctl_CMakeLists.txt
for file in ${path}roundctl/*
do
    cp $file roundctl/
done
mv roundctl_CMakeLists.txt roundctl/CMakeLists.txt

sed -i 's/call coll_getMinGapID(minGapID)/!call coll_getMinGapID(minGapID)/g' collimation.f90
sed -i 's/real(kind=fPrec), private, save :: emitnx0_dist    = zero/real(kind=fPrec), public, save :: emitnx0_dist    = zero/g' collimation.f90
sed -i 's/real(kind=fPrec), private, save :: emitny0_dist    = zero/real(kind=fPrec), public, save :: emitny0_dist    = zero/g' collimation.f90
sed -i 's/real(kind=fPrec), private, save :: emitnx0_collgap = zero/real(kind=fPrec), public, save :: emitnx0_collgap = zero/g' collimation.f90
sed -i 's/real(kind=fPrec), private, save :: emitny0_collgap = zero/real(kind=fPrec), public, save :: emitny0_collgap = zero/g' collimation.f90
sed -i 's/integer,          allocatable, private, save :: part_hit_pos(:)   /integer,          allocatable, public, save :: part_hit_pos(:)   /g' collimation.f90
sed -i 's/call shuffleLostParticles/!call shuffleLostParticles/g' collimation.f90
sed -i 's/if(ithick == 1) call synuthck/!if(ithick == 1) call synuthck/g' mod_particles.f90
