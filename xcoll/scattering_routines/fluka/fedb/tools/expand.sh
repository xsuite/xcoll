#!/bin/bash
#
# Vittorio Boccone -  22/06/10
#  boccone@cern.ch
#
# Modified by Roberto Versaci  22-09-2010
# Modified by Vittorio Boccone 02-03-2011 - add output file
#
# Fluka Project Builder 
#  replace the different includes
#  with the hardcoded text
#
inputfile=$1
#
base=${1%%.*} 
ext=${1#*.}
#
# if 2 arguments are passed use the second as destination
if [ "$#" = "2" ]  ; then
    outputfile=$2
else
    outputfile=$base"_exp."$ext
fi
#
if [ -n "${inputfile}" ] ; then
    if [ -e ${inputfile} ] ; then 
	echo -e "$(tput bold setaf 1)expand$(tput sgr0): Building the input file $(tput setaf 2)${outputfile}$(tput sgr0)"
	cp ${inputfile} ${outputfile}
	for name in $(awk -F'^#include' '{print $2}' ${inputfile}) ; do
	  newname=`echo "$name" | sed -e 's/\./\\\./g' | sed -e 's:\/:\\\/:g'`
	  echo -e "         Including $(tput setaf 4)${name}$(tput sgr0)"
	  if [ -e "$name" ] ; then
	      sed -e "/$newname/r $name" -e "/$newname/d" $outputfile > tmp.inp
	      mv tmp.inp $outputfile
	  else 
	      echo "File " $name " does not exist"
	  fi
	done
    else
	echo "File to be expanded does not exist"
    fi
else
    echo "No input file specified"
fi