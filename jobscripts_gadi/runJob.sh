#!/bin/bash
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=00:10:00
#PBS -l wd
#PBS -l mem=48GB

#Script: runJob.sh
#Purpose: sets up and runs a job on a gadi node. Sets relevant PBS environment variables
#         input file and relevant directory structure assumed to be set up in 
#         a script calling this one


module load intel-mkl/2021.2.0

export OMP_NUM_THREADS=48
export OMP_STACKSIZE=1000m

/scratch/d35/rh5686/H3+WorkDir/H3Plus

tracejob $PBS_JOBID > "$PBS_JOBID".trace

