#!/bin/bash

#SBATCH --job-name=form
#SBATCH --output=form-out-%j
#SBATCH --partition=research.q
#SBATCH --cpus-per-task=12
#SBATCH --array=0-X  # Definition of the number of arrays (X = n - 1)


# Start time.
echo "Start at $(date)".


# Storage of the combinations of names collected in "names_dir.tsv".

cd ../subfiles
readarray -t dir_array < <(awk 'NR>1 {print $1}' names_dir.tsv)

cd ../subscripts


# Execution of the subscripts according to the cell population selected and the SLURM array task ID.

Rscript U02_A_filter_names.R $(cd ../.. && pwd) ${CELL} ${dir_array[$SLURM_ARRAY_TASK_ID]}
wait

bash U02_B_create_form_list.sh $(cd ~/${CELL} && pwd) ${dir_array[$SLURM_ARRAY_TASK_ID]}
wait

mv ./*-out-* out_files


# End time.
echo "End at $(date)".


