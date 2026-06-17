#!/bin/bash

#SBATCH --job-name=comb
#SBATCH --output=comb-out-%j
#SBATCH --partition=research.q
#SBATCH --cpus-per-task=12


# Start time.
echo "Start at $(date)".


# Execution of the corresponding subscripts according to the selected cell type.

bash U01_A_rename_bams.sh
wait

Rscript U01_B_names_comb.R $(cd ../.. && pwd) ${CELL}

mv ./*-out-* out_files


# End time.
echo "End at $(date)".

