#!/bin/bash
#
#SBATCH --job-name=GridRun
#SBATCH --output=matlab_grid.“%j”.out
#SBATCH --error=matlab_grid.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < distMesh.m
