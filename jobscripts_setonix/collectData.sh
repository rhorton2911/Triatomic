#!/bin/bash


cd $MYSCRATCH/H3+WorkDir

readarray -t Rvals1 < R1vals.txt
readarray -t Rvals2 < R2vals.txt

for R1 in ${Rvals1[@]}
do 
   for R2 in ${Rvals2[@]}
   do
			islessthan=$(echo "${R1} < 2.0*${R2}" | bc -l)
			if (("$islessthan")); then
         R1Val=$(echo "${R1}" | awk '{printf "%.2f", $0}')
         R2Val=$(echo "${R2}" | awk '{printf "%.2f", $0}')

         cd R1${R1Val}R2${R2Val}
				 echo "$R1 $R2"
         for i in $(seq 1 4) 
         do 
            echo ${R1Val} ${R2Val} $(awk '$1=='${i} energies.txt) >> ../PEC.${i}
         done
				 cd ../
			fi
	 done
	 #Add a blank line between energies for different values of R1 to 
	 #allow for use of gnuplot
   for i in $(seq 1 4) 
   do 
      echo -e ' \n' >> ../PEC.${i}
   done
done
