#!/bin/bash


#=========================================
#	1) Definition of arguments
#=========================================

dir_bam=${1}
wait


#====================================================================================================
#	2) Displacement to the specified directory and storage of the whole list of BAM files.
#====================================================================================================

cd ${dir_bam}
list_f=( *.bam )


#================================================================================================================
#	3) Capture of the short identifier of each BAM file (second field separated by "_") in a for loop.
#================================================================================================================

for i in ${list_f[@]}
do
	ifin=$(echo "${i}" | cut -d'_' -f2)
	echo "${ifin}"

	# Change of the names of the files to the shorter version of the identifier.
	mv *${ifin}*.bam ${ifin}.bam
	mv *${ifin}*.bam.bai ${ifin}.bam.bai
	mv *${ifin}*.bam.txt ${ifin}.bam.txt

done


