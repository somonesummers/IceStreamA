#!/bin/bash
#
#SBATCH --job-name=TauRun
#SBATCH --output=matlab_test.“%j”.out
#SBATCH --error=matlab_test.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < tau_sandbox.m
