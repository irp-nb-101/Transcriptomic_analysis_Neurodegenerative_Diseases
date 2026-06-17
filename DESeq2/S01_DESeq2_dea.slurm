#!/bin/bash

#SBATCH --job-name=rmats
#SBATCH --output=rmats-out-%j
#SBATCH --partition=research.q
#SBATCH --mem=32G
#SBATCH --cpus-per-task=12


# Start time.
echo "Start at $(date)".

# Arguments definition.

# arg 1 --> cell_pop (cell population)
# arg 2 --> $PWD/ (current directory)
# arg 3 --> sample_type
# arg 4 --> pattern


# DESeq2 DEA script execution.

Rscript U01_A_RScript_DESeq_AD.R # {1} {2} {3} {4}


# Start time.
echo "Start at $(date)".









