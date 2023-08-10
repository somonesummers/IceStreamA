#!/bin/bash
#
#SBATCH --job-name=trimFiles
#SBATCH --output=matlab_reduce.“%j”.out
#SBATCH --error=matlab_reduce.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < trimData.m
