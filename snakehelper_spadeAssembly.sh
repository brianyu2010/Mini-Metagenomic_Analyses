#!/bin/bash

############################################
# Assembly helper via spades
############################################

scratch=$1
paired1=$2
paired2=$3
single1=$4
single2=$5
contigs=$6
scaffolds=$7
tool_dir=$8
tot_mem=$9	# this is in MB, need to get rid of right 3 zeros
threads=${10}
spades_output_dir=${11}

echo $scratch
cd $scratch
date
ls

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
mem=$( echo $tot_mem | rev | cut -c 4- | rev )
echo "memory used in Gb is "$mem

echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( wc -l < $paired1) / 4 ))

if [ -e $single1 ] && [ -e $single2 ]
then
  # zcat handles zip files directly, cat cannot
  cat $single1 $single2 > single.fastq
  echo -e 'Starting Assembly: Number of single end reads is: \t'$(( $( wc -l < single.fastq ) / 4 ))
else
  echo "No single end fastq reads detected"
fi

# Assembly -m stands for memory allowed
# Check if there are single end read files.
# Uneven coverage means that you start with small kmers and check mean read length 
# ~98 for nextseq samples and 140 for Hiseq samples
if [ -e $single1 ] && [ -e $single2 ]
then
  /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/spades.py" -k 21,33,55,77 -t $threads --careful --sc -m $mem \
    -1 $paired1 -2 $paired2 -s single.fastq -o $spades_output_dir
else
  /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/spades.py" -k 21,22,55,77 -t $threads --careful --sc -m $mem \
    -1 $paired1 -2 $paired2 -o $spades_output_dir
fi

cp $spades_output_dir/contigs.fasta $contigs
cp $spades_output_dir/scaffolds.fasta $scaffolds
mv $spades_output_dir/corrected/ corrected/

ls
echo
ls corrected/

mv corrected/*1.*.00.0_0.cor.fastq.gz correctedreads1.fastq.gz
mv corrected/*2.*.00.0_0.cor.fastq.gz correctedreads2.fastq.gz
mv corrected/*unpaired.00.0_0.cor.fastq.gz correctedreads_unpaired.fastq.gz

date
echo "Assembly Completed"

exit


