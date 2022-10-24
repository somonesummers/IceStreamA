#!/bin/bash
#
#SBATCH --job-name=Isa501
#SBATCH --output=matlab_Main.“%j”.out
#SBATCH --error=matlab_Main.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < MainHelper.m
