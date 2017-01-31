#!/bin/bash
# counts number of reads from each sample
# Usage: must be invoked in the folder with all the samples
# 
# Revision History:
# 2014.08.19 Created
# 2014.12.15 This is a MiSeq run so the folder names are different again

folder_list=($(ls))
for f in "${folder_list[@]}"
do
  cd $f
  echo $f
	# This is for Nextseq Names
	# sample_id=$(echo $f | cut -d_ -f 2,3)
	# This is for the new MiSeq folder names
	sample_id=$(echo $f | cut -d- -f 2,3)
  reads=$( zcat *R1_001.fastq.gz | wc -l )
  cd ..
  echo -e $sample_id'\t'$(($reads/4)) >> raw_reads.cnt
done

mv raw_reads.cnt ../results/

exit

