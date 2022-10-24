#!/bin/bash
#
#SBATCH --job-name=IsaOpt
#SBATCH --output=matlab_Opt.“%j”.out
#SBATCH --error=matlab_Opt.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < OptimizationRunner.m
