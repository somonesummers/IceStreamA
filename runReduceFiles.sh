#!/bin/bash
#
#SBATCH --job-name=trimFiles
#SBATCH --output=matlab_test.“%j”.out
#SBATCH --error=matlab_test.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < trimData.m
