#!/bin/bash
#SBATCH --job-name=hoomd-perf
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=RM-shared
#SBATCH --time=00:30:00
#SBATCH --account=see220002p

#test

N=500
singularity exec /ocean/projects/see220002p/shared/icomse_cpu.sif python perf.py $SLURM_JOB_ID $N $SLURM_NTASKS

