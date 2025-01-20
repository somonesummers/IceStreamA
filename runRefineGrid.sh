#!/bin/bash
#
#SBATCH --job-name=IsaMesh
#SBATCH --output=matlab_distmesh.“%j”.out
#SBATCH --error=matlab_distmesh.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < distMesh2.m
