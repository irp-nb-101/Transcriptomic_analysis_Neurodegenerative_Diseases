#!/bin/bash

#SBATCH --job-name=rmats
#SBATCH --output=rmats-out-%j
#SBATCH --partition=research.q
#SBATCH --mem=32G
#SBATCH --cpus-per-task=12
#SBATCH --array=0-X  # Definition of the number of arrays (X = n - 1)


# Start time.
echo "Start at $(date)".


# Arguments definition.

# arg 1 --> $PWD/ (current directory)
# arg 2 --> sample_type


# Creation of the array with the different DESeq2 results directories.

cd results/

deseq2_dir=( DES* )
current_dir_pre=${deseq2_dir[$SLURM_ARRAY_TASK_ID]}
current_dir=${current_dir_pre#DESeq2_}
echo "${current_dir}"

cd ..


# GSEA-GO and Volcano plot scripts execution.

Rscript U02_A_GSEA_GO_enrich.R $PWD/ ${current_dir}
wait

Rscript U02_B_Volcano_plots.R $PWD/ ${current_dir}


# Start time.
echo "Start at $(date)".




