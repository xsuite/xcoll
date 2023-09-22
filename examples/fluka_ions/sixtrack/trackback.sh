#!/bin/bash

LOSTIDS=($(awk 'NR>1{print $6}' aperture_losses.dat)) ###get the IDs of particles lost in the aperture file

echo " # PartID  TurnCreated  CollID  ParentID Proc Proc_turn1" >> ParticleSummary  ###create summary file header

TCPC=13 #collimator ID of the crystal

for i in "${LOSTIDS[@]}"
do
  turn=$(awk -v var="$i" '$2 == var {x=NR+3000} (NR<=x) {print}' fluka.log | grep -m1 '# FlukaIO: turn =' | awk '{print $5}')  ###look for the turn number (fifth word) in the 3000 lines after the appearence of the lost ID in fluka.log
  collID=$(awk -v var="$i" '$2 == var {x=NR+3000} (NR<=x) {print}' fluka.log | grep -m1 '# FlukaIO: turn =' | awk '{print $8}')  ###look for the collimator ID of where the particle was first created (eighth word) in the 3000 lines after the appearence of the lost ID in fluka.log
  parentID=$(awk -v var="$i" '$2 == var {print; printf "\n"}' fluka.log | head -n1 | awk '{print $3}') ###look for the parent ID (third word) in the first appearence of the lost ID in fluka.log
  crystalIDs=($(awk -v var=$TCPC '$1 == var {print; printf "\n"}' *fort.50 | awk -v var="$turn" '$21 == var {print; printf "\n"}' | awk -v var="$parentID" '$3 == var {print; printf "\n"}' | awk '{print $20}')) #list of processes with the right parentID, turn
  proc1=$(awk -v var=$TCPC '$1 == var {print; printf "\n"}' *fort.50 | awk '$21 == 1 {print; printf "\n"}' | awk -v var="$parentID" '$3 == var {print; printf "\n"}' | awk '{print $20}') #process in the first turn

  if [[ "${#crystalIDs[@]}" -gt 1 ]]; then
     proc="UNKNOWN"
  elif [[ "${#crystalIDs[@]}" -eq 0 ]]; then
     proc="INI"
  else
     proc=${crystalIDs[1]}
  fi
  
  if [ -z $proc1 ]; then
     proc1="INI"
  fi

  echo "$i  $turn  $collID  $parentID  $proc  $proc1" >> ParticleSummary

  awk -v var="$i" '$8 == var' lhcFTion001_allImpacts.dat >> ParticleImpacts_${i} ### writes down all the impacts of the particle
  awk -v var="$i" '$2 == var {print; printf "\n"}' fluka.log >> ParticleTrack_${i} ### writes down all the appearences of the lost particle in fluka.log
done
