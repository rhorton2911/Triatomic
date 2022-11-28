#!/bin/bash

#Script: runIso
#Purpose: runs the code over a grid of parameters describing the isosceles
#         H3++ molecule in a range of configurations.
#Running: uses absolute paths. Script can be run from anywhere on setonix. 
#         Assumes variables MYSOFWARE and MYSCRATCH have been defined 
#         to contain absolute paths of the user's /software and /scratch areas

if ! [[ -d /scratch/d35/rhorton/isosceleseTest/ ]]; then
   mkdir /scratch/d35/rhorton/isosceleseTest/
fi
cd $MYSOFTWARE/Triatomic/
cp ./H3Plus /scratch/d35/rhorton/isosceleseTest/

Rvals1=$(seq 0.2 0.4 7.0)
Rvals2=$(seq 0.2 0.4 7.0)

cd /scratch/d35/rhorton/isosceleseTest/
rm R1vals.txt
rm R2vals.txt
#Write energies to files for use in other scripts, avoids
#energy arrays being defined in mutiple places
for R1 in ${Rvals1[@]}
do
   echo "$R1" >> R1vals.txt
done
for R2 in ${Rvals2[@]}
do
   echo "$R2" >> R2vals.txt
done

cd $MYSOFTWARE/Triatomic/

#Declare array to save job IDs
declare -a idArray
idInd=0

for R1 in ${Rvals1[@]}
do 
	 #Value of R1 used in runJob.script
	 echo "$R1" > currR1Val
   jobId=$(sbatch --partition=work -J H3++R1=${R1} --error=H3++R1=${R1}.err $MYSOFTWARE/Triatomic/jobscripts/packjobscripts/runJobPack.script | cut -f 4 -d ' ')
   idArray[${idInd}]="$jobId"
   idInd=$((idInd+1))
	 rm currR1Val
done

#Setup list of depencencies for final cleanup script, slurm format is jobId:jobId:jobId:etc for all required jobs
rm /scratch/d35/rhorton/isosceleseTest/jobIdFile
dependencies=""
for id in "${idArray[@]}"
do
   echo "$id" >> /scratch/d35/rhorton/isosceleseTest/jobIdFile
   dependencies+=":${id}"
done

#The following two scripts either collect all the data produced into usable output files, or produce a list of failed jobs, depending on whether or
#not all jobs were successful.
sbatch --dependency=afterok"${dependencies}"  --partition=work -J H3++Cleanup $MYSOFTWARE/Triatomic/jobscripts/runDataCollect.script
sbatch --dependency=afterany"${dependencies}" --partition=work --time=00:01:00 -J H3++Errors $MYSOFTWARE/Triatomic/jobscripts/checkerrors.sh


