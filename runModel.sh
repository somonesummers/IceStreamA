#!/bin/bash
#
#SBATCH --job-name=IceStrA
#SBATCH --output=matlab_Main.“%j”.out
#SBATCH --error=matlab_Main.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < MainHelper.m
