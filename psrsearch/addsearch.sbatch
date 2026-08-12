#!/usr/bin/bash

#SBATCH --job-name=addsearch
#SBATCH --output=./logs/%x_%A_%a.out
#SBATCH --error=./logs/%x_%A_%a.out
#SBATCH --account=pawsey0272
#SBATCH --partition=work
#SBATCH --mem=16G
#SBATCH --time=6:00:00


export LD_LIBRARY_PATH=/software/projects/pawsey0272/zwan4817/lib:$LD_LIBRARY_PATH

addsearch /scratch/pawsey0272/zwan4817/J1646/20251225/add/links/*.sf -E add
