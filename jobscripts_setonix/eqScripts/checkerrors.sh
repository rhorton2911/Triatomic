#!/bin/bash

cd /scratch/d35/rhorton/H3+WorkDir

declare -a idArray
readarray -t idArray < jobIdFile

rm failedJobsFile 
for id in "${idArray[@]}"
do
	 exitCode=$(sacct -j "$id" | head -3 | tail -1 | xargs | rev | cut -c1-1 ) 
	 echo "$id $exitCode"
	 isNotZero=$(echo "$exitCode != 0" | bc)
	 if (("$isNotZero")); then
	    echo "$id $exitCode" >> failedJobsFile
	 fi
done


