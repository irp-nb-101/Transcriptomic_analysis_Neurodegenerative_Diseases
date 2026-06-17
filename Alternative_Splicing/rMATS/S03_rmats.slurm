#!/bin/bash

#SBATCH --job-name=rmats
#SBATCH --output=rmats-out-%j
#SBATCH --partition=research.q
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --array=0-X  # Definition of the number of arrays (X = n - 1)


# Start time.
echo "Start at $(date)".


# Activation of rMATS conda environment (modification of the source).

source </path/conda.sh>
conda activate rmats


# Definition of HOME directory and retrieval of the RefSeq mice (mm39) annotation data from NCBI (in GTF format). 

home_dir=$(cd ~ && pwd)
annotation_file=${home_dir}/Initial_data/proof_mm39.ncbiRefSeq.gtf


# Creation of all the possible combinations between each eleent of the control array and the corresponding test comparison.

cd control_bam_lists/
control_array=( * )

cd ../


full_array=( )
for i in ${control_array[@]}
do
	pre_file=${i%.txt}
	fin_file=$(echo "${pre_file}" | cut -d'_' -f3-)
	echo "${fin_file}"

	cd test_bam_lists/

	test_array=( *${fin_file}* )

	for j in ${test_array[@]}
	do
		concat=$(echo "${i}+${j}")
		echo "${concat}"

	full_array+=( ${concat} )
	
	done

	cd ../

done

echo "${full_array[@]}"

pre_ctrl_bam=$(echo "${full_array[$SLURM_ARRAY_TASK_ID]}" | cut -d'+' -f1)
pre_test_bam=$(echo "${full_array[$SLURM_ARRAY_TASK_ID]}" | cut -d'+' -f2)

echo "${pre_ctrl_bam}"
echo "${pre_test_bam}"

pre_sample=${pre_test_bam%.txt}
sample=${pre_sample#form_}

echo "${sample}"


# Selection of the control-test pair of BAM files.

ctrl_bam=$PWD/control_bam_lists/${pre_ctrl_bam}
test_bam=$PWD/test_bam_lists/${pre_test_bam}


# Definition of additional parameters for the rMATS execution.

sing_or_pair=${SP}
read_length=${L}



# Execution of rMATS with the selected parameters.

rmats.py --b1 ${ctrl_bam} --b2 ${test_bam} \
         --gtf ${annotation_file} \
	 --od $PWD/Ctrl-vs-${sample} \
	 --tmp $PWD/Ctrl-vs-${sample}/tmp \
         -t ${sing_or_pair} --readLength ${read_length} \
	 --variable-read-length \
         --nthread 12 &> rmats.log

wait


# End time.
echo "End at $(date)".


