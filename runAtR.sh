#!/bin/bash

cd diatomicCase
for R in 0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0
do 
   if ! [[ -d R=${R} ]]; then
      mkdir R=${R}
   fi

   #bc returns integers after division without "scale" argument
   RNuc=$(echo "scale=2; ${R}/2.0" | bc | awk '{printf "%.2f", $0}')

   sed s/RRR/${RNuc}/g input > R=${R}/input

   cd R=${R}
      export OMP_NUM_THREADS=10
      echo $OMP_NUM_THREADS
      #Run main executable stored in higher directory
      echo "Running: R=${R}"
      ../../H3Plus
   for i in $(seq 1 4) 
   do 
      echo ${R} $(awk '$1=='${i} energies.txt) >> ../PEC.${i}
   done
   cd ../ 
done

cd ../


