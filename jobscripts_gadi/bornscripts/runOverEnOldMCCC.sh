#!/bin/bash

#Script: runOverEn.sh

EVals1=$(seq 30.00 5.00 200.0)

scratchDir="/scratch/d35/rh5686/BornTestDirOldMCCC/"
#codeDir="~lhs573/code/MCCC-FN/build.gadi.intel/"
dataDir="/scratch/d35/rh5686/HeH+Code/"
if [ ! -d "$scratchDir" ]
then
   echo "Making directory: $scratchDir"
   mkdir "$scratchDir"
fi
#cd "$codeDir"
#cp main "$scratchDir"/
cd "$dataDir"
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

	 rm *ICS*
   #Remove previous pbs output files
   rm *.e*
   rm *.o*
	 rm *.trace

   l1=$(grep -n "ENERGY" data.in | cut -d ":" -f1) 

   #Replace lines with nuclear r coordinates
	 sed -i "${l1}s/:.*/: ${E} eV/" data.in  #Sed regex: .=any character, *=any number of occurences => .*= any number of any character (includes whitespace)

   #Submit PBS jobscript 
   job=$(qsub -N BornJobE"${E}" /g/data/d35/rh5686/Triatomic/jobscripts_gadi/bornscripts/runBornOldMCCC.sh)
   jobId=$( echo "$job" | cut -d "." -f1)
	 echo "$jobId" >> ../jobidlist

   idString+=":$jobId"
 
   cd ../
done

echo "$idString"
job=$(qsub -N BornCleanup -W depend=afterok"${idString}" /g/data/d35/rh5686/Triatomic/jobscripts_gadi/bornscripts/collectBornOldMCCC.sh)
jobId=$( echo "$job" | cut -d "." -f1)
echo "$jobId" > cleanupid
