#!/bin/bash

#Script: runEquilaterial.sh
#Purpose: runs H3+ structure code in equilaterial nuclear geometries
#         with varying internuclear distances, producing a H3+ potential
#         energy curve. Each value of R is a submitted as a separate job.


RVals1=$(seq 0.80 0.10 1.50)
RVals2=$(seq 1.55 0.05 1.75)
RVals3=$(seq 1.80 0.20 6.00)
RVals+=(${RVals1[@]} ${RVals2[@]} ${RVals3[@]})

scratchDir="/scratch/d35/rh5686/H3+WorkDir"
codeDir="/g/data/d35/rh5686/Triatomic"
if [ ! -d "$scratchDir" ]
then
   echo "Making directory: $scratchDir"
   mkdir "$scratchDir"
fi
cd "$codeDir"
cp H3Plus "$scratchDir"/
cp data.in "$scratchDir"/
cp input "$scratchDir"/
cd "$scratchDir"

rm RVals.txt
rm jobidlist
rm cleanupid
idString=""
for R in "${RVals[@]}"
do
   if [ ! -d "R${R}" ]
   then
      mkdir "R${R}"
   fi

   cp input "R${R}"/
   cp data.in "R${R}"/
   echo "${R}" >> RVals.txt
   cd "R${R}"

   #Remove previous pbs output files
   rm *.e*
   rm *.o*
	 rm *.trace

   #For an equilateral triangle, (r value of vertex)=(side length)/sqrt(3)
   r1=$( echo "scale=6; ${R}/sqrt(3)" | bc | awk '{printf "%.5f\n", $0}')  

   l1=$(grep -n "R1" input | cut -d ":" -f1) #linenumber of R1 variable
   l2=$(echo "$l1 + 1" | bc)
   l3=$(echo "$l2 + 1" | bc)

   #Replace lines with nuclear r coordinates
   sed -i "${l1}s/.*/${r1}   R1/" input
   sed -i "${l2}s/.*/${r1}   R2/" input
   sed -i "${l3}s/.*/${r1}   R3/" input 

   #Submit PBS jobscript 
   job=$(qsub -N H3+StructureR"${R}" /g/data/d35/rh5686/Triatomic/jobscripts_gadi/runJob.sh)
   jobId=$( echo "$job" | cut -d "." -f1)
	 echo "$jobID" >> ../jobidlist

   idString+=":$jobId"
 
   cd ../
done

echo "$idString"
job=$(qsub -N H3+Cleanup -W depend=afterok"${idString}" /g/data/d35/rh5686/Triatomic/jobscripts_gadi/collectData.sh)
jobId=$( echo "$job" | cut -d "." -f1)
echo "$jobId" > cleanupid
