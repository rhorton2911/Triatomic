#!/bin/bash

#Script: runEq
#Purpose: runs the code over a grid of parameters describing the equilateral
#         H3++ molecule with a range of internuclear distances
#Running: uses absolute paths. Script can be run from anywhere on setonix. 
#         Assumes variables MYSOFWARE and MYSCRATCH have been defined 
#         to contain absolute paths of the user's /software and /scratch areas

if ! [[ -d /scratch/d35/rhorton/equilateralCase/ ]]; then
   mkdir /scratch/d35/rhorton/equilateralCase
fi
cd $MYSOFTWARE/Triatomic/
cp ./H3Plus /scratch/d35/rhorton/equilateralCase/

Rvals=$(seq 0.2 0.4 7.0)

cd /scratch/d35/rhorton/equilateralCase/
rm Rvals.txt
for R in ${Rvals[@]}
do
   echo "$R" >> Rvals.txt
done
cd $MYSOFTWARE/Triatomic/

#Declare array to save job IDs
declare -a idArray
idInd=0

for R in ${Rvals[@]}
do
   RVal=$(echo "${R}" | awk '{printf "%.2f", $0}')
   if ! [[ -d /scratch/d35/rhorton/equilateralCase/R${RVal} ]]; then
      mkdir /scratch/d35/rhorton/equilateralCase/R${RVal}
   fi
  
   cat input > inputTemp 
   sed -i "s/R1/${RVal}/g" inputTemp
   sed -i "s/R2/${RVal}/g" inputTemp 
  
   cp data.in /scratch/d35/rhorton/equilateralCase/R${RVal}
   cp inputTemp /scratch/d35/rhorton/equilateralCase/R${RVal}/input
   rm inputTemp
  
   cd /scratch/d35/rhorton/equilateralCase/R${RVal}/
   #Run main executable stored in higher directory
   echo "Running: R=${RVal}"
  
   jobId=$(sbatch --partition=work -J H3++R=${R} --error=H3++R1${R1}R2${R2}.err $MYSOFTWARE/Triatomic/jobscripts/eqScripts/runJob.script | cut -f 4 -d ' ')
   idArray[${idInd}]="$jobId"
   idInd=$((idInd+1))
  
   cd /software/projects/d35/rhorton/Triatomic
done


#Setup list of depencencies for final cleanup script, slurm format is jobId:jobId:jobId:etc for all required jobs
rm /scratch/d35/rhorton/equilateralCase/jobIdFile
dependencies=""
for id in "${idArray[@]}"
do
   echo "$id" >> /scratch/d35/rhorton/equilateralCase/jobIdFile
   dependencies+=":${id}"
done

sbatch --dependency=afterok"${dependencies}"  --partition=work -J H3++Cleanup $MYSOFTWARE/Triatomic/jobscripts/eqScripts/runDataCollect.script
sbatch --dependency=afterany"${dependencies}" --partition=work --time=00:05:00 -J H3++Errors $MYSOFTWARE/Triatomic/jobscripts/eqScripts/checkerrors.sh


