#!/bin/bash
#
#SBATCH --job-name=IceStrA025
#SBATCH --output=matlab_Main.“%j”.out
#SBATCH --error=matlab_Main.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < MainHelper.m
