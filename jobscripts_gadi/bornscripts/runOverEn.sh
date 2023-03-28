#!/bin/bash

#Script: runOverEn.sh

EVals1=$(seq 30.00 5.00 200.0)

scratchDir="/scratch/d35/rh5686/BornTestDir/"
codeDir="/g/data/d35/rh5686/MCCC-FN-Triatomic/build.gadi.intel/"
if [ ! -d "$scratchDir" ]
then
   echo "Making directory: $scratchDir"
   mkdir "$scratchDir"
fi
cd "$codeDir"
cp main "$scratchDir"/
cd ../
cp data.in "$scratchDir"/
cd "$scratchDir"

rm EVals.txt
rm jobidlist
rm cleanupid
idString=""
for E in ${EVals1[@]}
do
   if [ ! -d "E${E}" ]
   then
      mkdir "E${E}"
   fi

   cp data.in "E${E}"/
   echo "${E}" >> EVals.txt
   cd "E${E}"

   #Remove previous pbs output files
   rm *.e*
   rm *.o*
	 rm *.trace

   l1=$(grep -n "ENERGY" data.in | cut -d ":" -f1) 

   #Replace lines with nuclear r coordinates
	 sed -i "${l1}s/:.*/: ${E} eV/" data.in  #Sed regex: .=any character, *=any number of occurences => .*= any number of any character (includes whitespace)

   #Submit PBS jobscript 
   job=$(qsub -N BornJobE"${E}" /g/data/d35/rh5686/Triatomic/jobscripts_gadi/bornscripts/runBorn.sh)
   jobId=$( echo "$job" | cut -d "." -f1)
	 echo "$jobId" >> ../jobidlist

   idString+=":$jobId"
 
   cd ../
done

echo "$idString"
job=$(qsub -N BornCleanup -W depend=afterok"${idString}" /g/data/d35/rh5686/Triatomic/jobscripts_gadi/bornscripts/collectBorn.sh)
jobId=$( echo "$job" | cut -d "." -f1)
echo "$jobId" > cleanupid
