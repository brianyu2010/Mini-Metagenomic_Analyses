#!/bin/bash

############################################
# Spades assembly to be run on bigmem nodes
# Or these can be run for bulk sequences
# 2016.09.23 This file is updated to use
#            metaspades ie. spades --meta so no --careful flag
# 2017.01.13 Changed input arguments from one single
#            argument to two. In the future all
#            should use YAML input file
############################################

# All input should be zipped files

scratch=$1
# pairedInterlaced=$2
paired1=$2
paired2=$3
# single1=$4
# single2=$5
code_dir=$4
tool_dir=$5    # for SPAdes it's /local10G/brianyu/tools/
tot_mem=$6     # this is in MB, need to get rid of right 3 zeros
threads=$7
spades_output_dir=$8

echo $scratch
# 2015.12.03 Not sure if you should cd into the scratch directory
# you cannot cd to scratch if you pass in input (not input_on_scratch)
# This is because once you cd to scratch, you cannot access ILxxxx..../P1....
# On the other hand, if you don't cd to scratch try not to create any new files.
cd $scratch
cp $code_dir/fastq_pairedEndToInterlace.py $scratch

date
pwd
ls
echo

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
mem=$( echo $tot_mem | rev | cut -c 4- | rev )
echo "memory used in Gb is "$mem

# making an interlased fastq file
# python fastq_pairedEndToInterlace.py $paired1 $paired2 pairedInterlaced.fastq

echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( cat $paired1 | wc -l ) / 4 ))

# Taking care of single reads; No need because metaspades cannot use single reads
# if [ -e $single1 ] && [ -e $single2 ]
# then
  # zcat handles zip files directly, cat cannot       
  # cat $single1 $single2 > single.fastq
  # echo -e 'Starting Assembly: Number of single end reads is: \t'$(( $( wc -l < single.fastq ) / 4 ))
# else
  # echo "No single end fastq reads detected"
# fi

# Checking if both single files are given
# if [ -e $single1 ] && ! [ -e $single2 ]
# then
  # echo "Only "$single1" found"
# fi

# if [ -e $single2 ] && ! [ -e $single1 ]
# then
  # echo "Only "$single2" found"
# fi

# Assembly -m stands for memory allowed
# Check if there are single end read files.
# Spades can accept .fastq.gz files
# uneven coverage means you start with smaller kmers and check your mean 
# read length so you are not using too high of kmers.
# Nextseq is around 98bp mean read length, hiseq is ~140bp
#
# Could also assert --only-assembler but for now did not
# must be one paired end input file in interlaced format. No single reads


if [ -d $spades_output_dir ]
then
  /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/metaspades.py" -k 21,33,55,77,91 -t $threads --continue -m $mem -1 $paired1 -2 $paired2 -o $spades_output_dir
else
  /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/metaspades.py" -k 21,33,55,77,91 -t $threads -m $mem -1 $paired1 -2 $paired2 -o $spades_output_dir
fi


ls | wc -l
date
echo "Assembly Completed"

exit


