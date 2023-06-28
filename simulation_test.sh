#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10000M   # memory per CPU core
#SBATCH --nodes=1 --ntasks=6 --mem=12G --gpus=1




# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load cuda
module load python/3.7

matlab -nodisplay -nojvm -nosplash -r "clear all;DDPM_microstructures_1chanel_64"