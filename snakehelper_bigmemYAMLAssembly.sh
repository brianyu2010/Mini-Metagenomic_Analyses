#!/bin/bash

############################################
# Spades assembly to be run on bigmem nodes
# 2017.01.12 Edited to use YAML file and combined
#            subsample contigs as --untrusted-contigs
############################################

# All input should be zipped files

scratch=$1
yaml_file=$2
tool_dir=$3
tot_mem=$4     # this is in MB, need to get rid of right 3 zeros
threads=$5
spades_output_dir=$6

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

# Assembly -m stands for memory allowed
# Spades can accept .fastq.gz files
# uneven coverage means you start with smaller kmers and check your mean 
# read length so you are not using too high of kmers.
# Nextseq is around 98bp mean read length, hiseq is ~140bp

# 2017.01.22 removed --careful due to time (performance) considerations

if [ -d $spades_output_dir ]
then
  /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/spades.py" --only-assembler -k 21,33,55,77,99 -t $threads --continue --sc -m $mem --dataset $yaml_file -o $spades_output_dir
else
  /usr/bin/time -v python $tool_dir"Spades-3.9.0/bin/spades.py" --only-assembler -k 21,33,55,77,99 -t $threads --sc -m $mem --dataset $yaml_file -o $spades_output_dir
fi

ls | wc -l
date
echo "Assembly Completed"

exit


