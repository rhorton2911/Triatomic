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

Rvals1=$(seq 0.2 0.2 7.0)
Rvals2=$(seq 0.2 0.2 7.0)

cd /scratch/d35/rhorton/isosceleseTest/
rm R1vals.txt
rm R2vals.txt
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
   for R2 in ${Rvals2[@]}
   do
			islessthan=$(echo "${R1} < 2.0*${R2}" | bc -l)
			if (("$islessthan")); then
         R1Val=$(echo "${R1}" | awk '{printf "%.2f", $0}')
         R2Val=$(echo "${R2}" | awk '{printf "%.2f", $0}')

         if ! [[ -d /scratch/d35/rhorton/isosceleseTest/R1${R1Val}R2${R2Val} ]]; then
            mkdir /scratch/d35/rhorton/isosceleseTest/R1${R1Val}R2${R2Val}
         fi

				 cat input > inputTemp 
         sed -i "s/R1/${R1Val}/g" inputTemp
         sed -i "s/R2/${R2Val}/g" inputTemp 

         cp data.in /scratch/d35/rhorton/isosceleseTest/R1${R1Val}R2${R2Val}/
         cp inputTemp /scratch/d35/rhorton/isosceleseTest/R1${R1Val}R2${R2Val}/input
				 rm inputTemp

         cd /scratch/d35/rhorton/isosceleseTest/R1${R1Val}R2${R2Val}/
         #Run main executable stored in higher directory
         echo "Running: R1=${R1Val} and R2=${R2Val}"

         jobId=$(sbatch --partition=work -J H3++R1${R1}R2${R2} --error=H3++R1${R1}R2${R2}.err $MYSOFTWARE/Triatomic/jobscripts/runJob.script | cut -f 4 -d ' ')
         idArray[${idInd}]="$jobId"
         idInd=$((idInd+1))

         cd /software/projects/d35/rhorton/Triatomic
			fi
	 done
done


#Setup list of depencencies for final cleanup script, slurm format is jobId:jobId:jobId:etc for all required jobs
rm /scratch/d35/rhorton/isosceleseTest/jobIdFile
dependencies=""
for id in "${idArray[@]}"
do
   echo "$id" >> /scratch/d35/rhorton/isosceleseTest/jobIdFile
   dependencies+=":${id}"
done

sbatch --dependency=afterok"${dependencies}"  --partition=work -J H3++Cleanup $MYSOFTWARE/Triatomic/jobscripts/runDataCollect.script
sbatch --dependency=afterany"${dependencies}" --partition=work --time=00:01:00 -J H3++Errors $MYSOFTWARE/Triatomic/jobscripts/checkerrors.sh


