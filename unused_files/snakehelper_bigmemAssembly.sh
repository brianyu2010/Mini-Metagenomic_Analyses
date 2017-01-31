#!/bin/bash

############################################
# Spades assembly to be run on bigmem nodes
############################################

# All input should be zipped files

scratch=$1
paired1=$2
paired2=$3
single=$4
tool_dir=$5
tot_mem=$6     # this is in MB, need to get rid of right 3 zeros
threads=$7
spades_output_dir=$8

echo $scratch
# 2015.12.03 Not sure if you should cd into the scratch directory
# you cannot cd to scratch if you pass in input (not input_on_scratch)
# This is because once you cd to scratch, you cannot access ILxxxx..../P1....
# On the other hand, if you don't cd to scratch try not to create any new files.
# cd $scratch
date
pwd
# ls
echo

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
mem=$( echo $tot_mem | rev | cut -c 4- | rev )
echo "memory used in Gb is "$mem

echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( zcat $paired1 | wc -l ) / 4 ))

if [ -e $single ]
then
  echo -e 'Starting Assembly: Number of single end reads is: \t'$(( $( zcat $single | wc -l ) / 4 ))
else
  echo "No single end fastq reads detected"
fi

# Assembly -m stands for memory allowed
# Check if there are single end read files.
# Spades can accept .fastq.gz files
# uneven coverage means you start with smaller kmers and check your mean 
# read length so you are not using too high of kmers.
# Nextseq is around 98bp mean read length, hiseq is ~140bp
if [ -e $single ]
then
  if [ -d $spades_output_dir ]
  then
    /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/spades.py" --only-assembler -k 21,33,55,77 -t $threads --careful --continue --sc -m $mem -1 $paired1 -2 $paired2 -s $single -o $spades_output_dir
  else
    /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/spades.py" --only-assembler -k 21,33,55,77 -t $threads --careful --sc -m $mem -1 $paired1 -2 $paired2 -s $single -o $spades_output_dir
  fi

else

  if [ -d $spades_output_dir ]
  then
    /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/spades.py" --only-assembler -k 21,33,55,77 -t $threads --careful --continue --sc -m $mem -1 $paired1 -2 $paired2 -o $spades_output_dir
  else
    /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/spades.py" --only-assembler -k 21,33,55,77 -t $threads --careful --sc -m $mem -1 $paired1 -2 $paired2 -o $spades_output_dir
  fi
fi

ls | wc -l
date
echo "Assembly Completed"

exit


