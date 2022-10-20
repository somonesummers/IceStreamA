#!/bin/bash
#
#SBATCH --job-name=IceStrA
#SBATCH --output=matlab_test.“%j”.out
#SBATCH --error=matlab_test.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < lib/gridSandbox.m
