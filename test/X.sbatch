#!/bin/bash
#SBATCH --job-name=ReO3press
#SBATCH --mail-user=lianyl
#SBATCH --partition=NEU
#SBATCH --nodes=1
#SBATCH --time=40:00:00
#SBATCH --array=1-6
#SBATCH --output=ReO3_press_%A_%a.log
#SBATCH --error=ReO3_press_%A_%a.err

echo "--------------- PROGRAM OUTPUT ----------------"

/csns_workspace/CSG/lianyl/softwares/julia/bin/julia   test_ReO3_pressure.jl

echo "--------------- PROGRAM COMPLETED--------------"

exit
