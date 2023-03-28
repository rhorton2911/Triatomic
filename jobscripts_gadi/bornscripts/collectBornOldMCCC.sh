#!/bin/bash
#PBS -q normal
#PBS -l ncpus=1
#PBS -l walltime=00:00:30
#PBS -l wd

#Script: collectBorn
#Purpose: collects two electron energies computed for equilateral geometries
#         and produces a file of potential energy curves.


cd /scratch/d35/rh5686/BornTestDirOldMCCC
rm BornICS.txt

readarray -t EVals < EVals.txt

echo "E(a.u)    Nf (1 -> N_max)"
for E in "${EVals[@]}"
do
   cd E"${E}"
   echo "E= ${E}"

	 #Collect born cross section
	 fileen=$(printf "%1.4E" ${E})
	 filename="ICS_"
	 filename+="${fileen}"
   line=$(grep -n "TICS" "${filename}" | cut -d ":" -f1)
   startval=$(echo "${line} + 1" | bc)
	 fullline=$(sed "${startval}q;d" "${filename}" )

   string="${E}       ${fullline}"   
	 echo "${string}" >> ../BornICS.txt

   cd ../
done
