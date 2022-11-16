#!/bin/bash
#SBATCH --job-name=hoomd
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=RM-shared
#SBATCH --time=00:30:00
#SBATCH --account=see220002p

T=1.4
singularity exec /ocean/projects/see220002p/shared/icomse_cpu.sif python test-cpu.py $SLURM_JOB_ID $T

