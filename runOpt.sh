#!/bin/bash
#
#SBATCH --job-name=OptThm025
#SBATCH --output=matlab_Opt.“%j”.out
#SBATCH --error=matlab_Opt.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < OptimizationRunner.m
