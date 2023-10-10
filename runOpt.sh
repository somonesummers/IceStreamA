#!/bin/bash
#
#SBATCH --job-name=Opt035
#SBATCH --output=matlab_Opt.“%j”.out
#SBATCH --error=matlab_Opt.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < OptimizationRunner.m
