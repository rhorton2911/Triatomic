#!/bin/bash

ls PEC.* &>/dev/null && rm PEC.* #Remove files from previous run
cp ./assignment_2_files/PEC.* ./

for R in 0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0
do 
	 if ! [[ -d R=${R} ]]; then
      mkdir R=${R}
	 fi
   sed s/RRR/${R}/ input > R=${R}/input

   cd R=${R}
	 #Run main executable stored in above directory
	 ../H2Plus
   for i in $(seq 1 4) 
   do 
      echo ${R} $(awk '$1=='${i} energies.txt) >> ../PEC.${i}
   done
   cd ../ 
done


