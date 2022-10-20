#!/bin/bash
#
#SBATCH --job-name=IceStrA
#SBATCH --output=matlab_test.“%j”.out
#SBATCH --error=matlab_test.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < MainHelper.m
