#!/bin/bash


cd $MYSCRATCH/equilateralCase

readarray -t Rvals < Rvals.txt

for R in ${Rvals[@]}
do 
   RVal=$(echo "${R}" | awk '{printf "%.2f", $0}')
 
   cd R${RVal}
 	 echo "$R"
   for i in $(seq 1 4) 
   do 
      echo ${RVal} $(awk '$1=='${i} energies.txt) >> ../PEC.${i}
   done
 	 cd ../
 
	 #Add a blank line between energies for different values of R1 to 
	 #allow for use of gnuplot
   for i in $(seq 1 4) 
   do 
      echo -e ' \n' >> ../PEC.${i}
   done
done
