#!/bin/bash
#
#SBATCH --job-name=1Rise02DhDt
#SBATCH --output=matlab_Main.“%j”.out
#SBATCH --error=matlab_Main.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=52:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < MainHelper.m
