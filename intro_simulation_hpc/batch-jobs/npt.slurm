#!/bin/bash
#SBATCH --job-name=hoomd
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --partition=GPU-shared
#SBATCH --gres=gpu:v100-16:1
#SBATCH --time=00:30:00
#SBATCH --account=see220002p

echo "testing lj-npt on one gpu"
T=1.3
singularity exec --nv /ocean/projects/see220002p/shared/icomse_gpu.sif python lj-npt.py $SLURM_JOB_ID $T

