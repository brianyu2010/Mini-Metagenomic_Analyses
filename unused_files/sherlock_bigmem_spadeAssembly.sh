#!/bin/bash

#SBATCH --job-name=Assembly  # job name   
#SBATCH --qos=bigmem
##SBATCH --qos=long   
#SBATCH --ntasks=1                      # number of nodes to allocate per job  
#SBATCH --cpus-per-task=5               # cpus (threads) per task
#SBATCH --partition=bigmem              # partition (queue) to use.
##SBATCH --partition=normal
#SBATCH --mem=240000                    # total memory in MB  
#SBATCH --time=16:00:00                 # wall clock time (d-H:MM:SS or minutes) 

## This is a comment and not recognized by the cluster  

## script assumes length of seedfile provided by command line with "--array=1-<file length>" 

# File Description:
# This file is run on brianyu@sherlock.stanford.edu 
# scratch space is on $SCRATCH/brianyu/  
# PI space is on $PI_HOME/brianyu/

# Revision History
# 2015.07.12 Created 

############################################
# Assembly helper via spades
############################################

# All inputs are zipped files

paired1=$1
paired2=$2
single=$3
spades_output_dir=$4
tot_mem=240000	# this is in MB, need to get rid of right 3 zeros
threads=5

echo $LOCAL_SCRATCH
cd $LOCAL_SCRATCH
date
ls

module load python/2.7.5/
tool_dir=$HOME/tools/

# cp $paired1 p1.fastq.gz
# cp $paired2 p2.fastq.gz

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
if [ -e $single ]
then
  /usr/bin/time -v python $tool_dir"Spades-3.5.0/bin/spades.py" --only-assembler -k 33,55,77,99 -t $threads --careful --continue --sc -m $mem -1 $paired1 -2 $paired2 -s $single -o $spades_output_dir
else
  /usr/bin/time -v python $tool_dir"Spades-3.5.0/bin/spades.py" --only-assembler -k 33,55,77,99 -t $threads --careful --sc -m $mem -1 $paired1 -2 $paired2 -o $spades_output_dir
fi

# Copy files back

ls
# spades should create output directly in $SCRATCH
# cp -r $spades_output_dir/ $PI_HOME/brianyu/results/

# gzip $paired1
# gzip $paired2

date
echo "Assembly Completed"

exit


