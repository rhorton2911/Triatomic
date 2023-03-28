#!/bin/bash
#PBS -q normal
#PBS -l ncpus=1
#PBS -l walltime=00:05:00
#PBS -l wd
#PBS -l mem=2GB
#PBS -P iy23
#PBS -lstorage=scratch/d35+gdata/d35

#Script: runJob.sh
#Purpose: sets up and runs a job on a gadi node. Sets relevant PBS environment variables
#         input file and relevant directory structure assumed to be set up in 
#         a script calling this one


module load intel-mkl/2021.2.0
module load openmpi/4.1.4

export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1000m

#mpirun -n 1 --map-by ppr:1:socket -bind-to socket /scratch/d35/rh5686/H3+WorkDir/H3Plus
mpirun -n 1 /scratch/d35/rh5686/BornTestDirOldMCCC/main

tracejob $PBS_JOBID > "$PBS_JOBID".trace

