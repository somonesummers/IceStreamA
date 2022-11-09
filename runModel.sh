#!/bin/bash
#
#SBATCH --job-name=1Rise02D
#SBATCH --output=matlab_Main.“%j”.out
#SBATCH --error=matlab_Main.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < MainHelper.m
