#!/bin/bash

if [ $# -ne 1 ]
then
  echo "This script needs exactly 1 argument: the path to the original SixTrack source."
  exit 1
fi

path=${1%/}/

cd FORTRAN_src

for f in core_tools constants strings mod_alloc common_modules string_tools mod_units mod_particles mod_meta mod_time bouncy_castle libcrlibm libroundctl coll_jawfit coll_common coll_db mod_ranlux mod_funlux coll_crystal coll_k2 coll_dist collimation
do
    cp ${path}${f}.f90 .
done

sed -i 's/integer, private, save :: c_ix  /integer, public,  save :: c_ix  /g' collimation.f90
sed -i 's/real(kind=fPrec), private, save :: emitnx0_dist    = zero/real(kind=fPrec), public, save :: emitnx0_dist    = zero/g' collimation.f90
sed -i 's/real(kind=fPrec), private, save :: emitny0_dist    = zero/real(kind=fPrec), public, save :: emitny0_dist    = zero/g' collimation.f90
sed -i 's/real(kind=fPrec), private, save :: emitnx0_collgap = zero/real(kind=fPrec), public, save :: emitnx0_collgap = zero/g' collimation.f90
sed -i 's/real(kind=fPrec), private, save :: emitny0_collgap = zero/real(kind=fPrec), public, save :: emitny0_collgap = zero/g' collimation.f90
sed -i 's/integer,          allocatable, private, save :: part_hit_pos(:)   /integer,          allocatable, public, save :: part_hit_pos(:)   /g' collimation.f90
sed -i 's/call shuffleLostParticles/!call shuffleLostParticles/g' collimation.f90
sed -i 's/if(ithick == 1) call synuthck/!if(ithick == 1) call synuthck/g' mod_particles.f90
