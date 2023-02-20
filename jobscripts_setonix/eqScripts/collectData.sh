#!/bin/bash


#cd $MYSCRATCH/H3+WorkDir
#
#readarray -t Rvals < Rvals.txt
#
#for R in ${Rvals[@]}
#do 
#   RVal=$(echo "${R}" | awk '{printf "%.2f", $0}')
# 
#   cd R${RVal}
# 	 echo "$R"
#   for i in $(seq 1 4) 
#   do 
#      echo ${RVal} $(awk '$1=='${i} energies.txt) >> ../PEC.${i}
#   done
# 	 cd ../
# 
#	 #Add a blank line between energies for different values of R1 to 
#	 #allow for use of gnuplot
#   for i in $(seq 1 4) 
#   do 
#      echo -e ' \n' >> ../PEC.${i}
#   done
#done

cd $MYSCRATCH/H3+WorkDir

readarray -t RVals < RVals.txt

mv PEC PEC_old
mv PEC_1e PEC_1e_old

echo "R(a.u)          E1(a.u)           E2(a.u)            E3(a.u)            E4(a.u)            E5(a.u)            E6(a.u)             E7(a.u)           E8(a.u)" >> PEC
echo "R(a.u)          E1(a.u)           E2(a.u)            E3(a.u)            E4(a.u)            E5(a.u)            E6(a.u)             E7(a.u)           E8(a.u)" >> PEC_1e
for R in "${RVals[@]}"
do
   cd R"${R}"
   echo "R= ${R}"

	 #Collect two electron structure results
   line=$(grep -n "E(a.u)" 2eenergies.txt | cut -d ":" -f1)
   firsteight=$(echo "${line} + 8" | bc)
   linenums=$(seq "${line}" 1 "${firsteight}")
   string="${R}   "

   startval=$(echo "${line} + 1" | bc)
   endval=firsteight
   for ((num=startval; num<=endval; num++)) 
   do
      en=$(sed -n "${num}p" 2eenergies.txt | awk '{print $NF}')
      string="${string}       ${en}"   
   done
   echo "${string}" >> ../PEC
 
	 #Collect one electron structure results
   line=$(grep -n "E(Ha)" 1eenergies.txt | cut -d ":" -f1)
   firsteight=$(echo "${line} + 8" | bc)
   linenums=$(seq "${line}" 1 "${firsteight}")
   string="${R}   "

   startval=$(echo "${line} + 1" | bc)
   endval=firsteight
   for ((num=startval; num<=endval; num++)) 
   do
      en=$(sed -n "${num}p" 1eenergies.txt | awk '{print $NF}')
      string="${string}       ${en}"   
   done
   echo "${string}" >> ../PEC_1e

   cd ../
done


