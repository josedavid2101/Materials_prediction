#!/bin/bash

#SBATCH --time=25:00:00   # walltime
#SBATCH --ntasks=5   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2000M   # memory per CPU core


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load matlab
module load gcc/12

matlab -nodisplay -nojvm -nosplash -r "clear all;main_get_odf_fourier_sc_1"