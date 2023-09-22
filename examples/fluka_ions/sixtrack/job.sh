#!/bin/bash

if [ $# -ne 1 ]
then
  echo "The script job.sh needs 1 argument: ./job.sh seed"
  exit 1
fi

seed=$1


cd ${path}/Results
study=TEST_s${seed}
if [ -d $study ] && [ -f ${study}/aperture_losses.dat ]
then
  echo "Results already exist for study "${study}"! Aborting..."
  exit 1
fi
mkdir $study
cd $study
ln -fns ${path}/Input/fort.2_FLUKA_${beam}${plane} fort.2
ln -fns ${path}/Input/fort.8_${beam} fort.8
ln -fns ${path}/Input/fort3.limi_FLUKA_${beam} fort3.limi
ln -fns ${path}/Input/insertion_${beam}.txt insertion.txt
ln -fns ${path}/Input/lhc_coupling_${beam}_exp.inp lhc_coupling_exp.inp
ln -fns ${path}/Input/relcol.dat_${beam} relcol.dat
ln -fns ${path}/Input/CollDB-RunII_${beam}.dat CollDB-RunII.dat

cp ${path}/Input/fort.3_FLUKA_${beam} fort.3
sed -i 's#%STUDY%#'${study}'   (Path: '${path}')#g' fort.3
thisseed=$(( $seed + 1574563 ))
sed -i 's/%SEED%/'${thisseed}'/g' fort.3

cp ${path}/Input/gpdist_input.txt_FLUKA_${beam}${plane} gpdist_input.txt
thisseed=$(( $seed + 38745 ))
sed -i 's/%SEED%/'${thisseed}'/g' gpdist_input.txt

# Executables
cp -P ${path}/Input/sixtrack .
cp -P ${path}/Input/flukaserver .
cp -P ${path}/Input/gpdist.exe .

# make a temp dir on the node and run SixTrack there
cd ${path}/Results
temp=$( mktemp -d )
cp -rP $study ${temp}/
cd ${temp}/$study

# Make initial distribution
./gpdist.exe gpdist_input.txt

# FLUKA environment
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/9.2.0/x86_64-centos7/setup.sh
export FLUSPC=/eos/project/f/flukafiles/fluka-coupling
export FLUKA=${FLUSPC}/fluka-4-3.3.Linux-gfor9
export CPLPATH=${FLUSPC}/fluka_coupling
export FEDBPATH=${FLUSPC}/fedb_coupling

# Declare network host
hostname > network.nfo

# Run the server in background:
${FLUKA}/bin/rfluka lhc_coupling_exp -e flukaserver -M 1 >> ${path}/Results/${study}/flukaserver.log 2>&1 &

# Wait until server is initialized:
while [ ! -f network.nfo ]; do
  echo "server still not running, waiting..."
  sleep 10
done

# wait until network.nfo declares both host and port
while [ `\wc -l network.nfo | \awk '{print ($1)}'` -ne 2 ] ; do
  echo "host/port not fully declared, waiting..."
  sleep 10
  # Check FLUKA is still running on the background. If not, close the job.
  if [ ` \ps -u $(whoami) | \grep rfluka | \wc -l ` -eq 0 ]; then
      echo "Something went wrong: FLUKA exited before opening the port!" 1>&2
      exit 1
  fi
done

# Run SixTrack
./sixtrack  > ${path}/Results/${study}/six.out 2>&1

sleep 20
echo " - Tracker done"

# Wait until server background process ends
wait
echo " - Server done"

# Retrieve results
echo "# turn collID interactiontype nImpacts" > fluka_absorptions_per_turn.dat
for f in lhc_coupling_exp*_lossMap.dat
do
  tail -n +2 $f | awk -- '{print $10" "$1" "$8;}' >> tempimpacts
done
sort -n tempimpacts | uniq -c | awk -- '{print $2" "$3" "$4" "$1;}' >> fluka_absorptions_per_turn.dat

cp aperture_losses.dat fluka_absorptions_per_turn.dat lhc_coupling_exp*_toucMap.dat ${path}/Results/${study}/
#rm CollPositions.dat MaterialInformation.txt align_error.dat collgaps.dat collsettings.dat colltrack.out file_units.log fluka.log fluka_isotope.log fort.{4,9,18,19,20,208,21,26,31} linopt_coupled.dat linopt_dump.dat orbitchecking.dat pencilbeam_distr.dat random_seeds.dat sigmasettings.out sim_meta.dat sim_time.dat singletrackfile.dat survival.dat tracks2.dat

