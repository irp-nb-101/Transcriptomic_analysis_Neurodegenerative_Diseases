#!/bin/bash


#=========================================
#	1) Definition of arguments
#=========================================

bams_path=${1}
dir=${2}

wait


#============================================================================================================================================
#	2) Specification of the lists with the names of the filtered BAM file names and the path where all of the BAM files are located
#============================================================================================================================================

list_test_ids=$(cat ../subfiles/filtered_bam_names/filt_${dir}.tsv)

if [[ ${dir} == *_f ]]; then
bams_dir=${bams_path}/BAM_pv_females

else
bams_dir=${bams_path}/BAM_pv_males

fi

out_file=form_${dir}.txt
> pre_file_${dir}.txt


#========================================================================================
#       3) Creation of the file with the complete BAM file paths separated by commas
#========================================================================================

for i in ${list_test_ids}
do
	if ! [ ${i} == "FileName1" ];
	then
		echo "${i}"
		echo "${bams_dir}/${i}.bam," >> pre_file_${dir}.txt
	fi
done

paste -sd '' pre_file_${dir}.txt | sed 's/.$//' > ${out_file}


if [[ ${out_file} == *Control* ]]; then
mv ${out_file} ../rMATS/control_bam_lists/
else
mv ${out_file} ../rMATS/test_bam_lists/
fi
rm pre_file_${dir}.txt

