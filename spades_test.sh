#!/bin/bash

#SBATCH --job-name=Assemble 	# job name
#SBATCH --ntasks=1           # number of nodes to allocate per job
#SBATCH --cpus-per-task=12		# cpus (threads) per task
#SBATCH --partition=bigmem		# partition (queue) to use.
##SBATCH --mem-per-cpu=5333		# memory request per cpu (threads) in MB
#SBATCH --mem=120000			# total memory
##SBATCH --time=72:00:00	# wall clock time (d-H:MM:SS or minutes)

## This is a comment and not recognized by the cluster

## script assumes length of seedfile provided by command line with "--array=1-<file length>"

# File Description:
# This file is run on /local10G/brianyu/ and is called by segmentReads4Assembly.sh

# Revision History
# 2014.05.16 Created
# 2014.10.27 Edited to use DNA clust first and then spade and do quast. Also keep kmer frequency
# 2014.12.11 Added ProcessClusteredFastq.py to re-order paired end reads
# 2014.12.12 Removed seedfile as an entry. Instead, experiment folder, similarity, kmer, and output_folder are passed in
#						 This file now must be called by another shell script which passes in these values
# 2014.12.14 Actually Added back the functionality of passing in a seed file
# 2014.12.28 Adapted for single cell cluster at /local10G/
# 2015.01.07 Changed to Version 2 where no reference is available and included trimming step.
#            Still does not include blast step after assembly
# 2015.04.13 Added blast function and removed Quast. This is for YNP environmental samples.

date
hostname


##################################
# User Set Directories
##################################

paired1=$1
paired2=$2
numreads=$3
output_dir=$4
#single1=$3
#single2=$4
#contigs=$6
#scaffolds=$7
#logfile=$8
#tool_dir=$9
#tot_mem=${10}	# this is in MB, need to get rid of right 3 zeros

scratch=$LOCAL_SATA

echo $scratch
cd $scratch
date

cp $paired1 paired1.fastq
cp $paired2 paired2.fastq
cp /local10G/brianyu/SnakeMake_Projects/DNA_Assembly/randomizeFastq.py .

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
#mem=$( echo $tot_mem | rev | cut -c 4- | rev )
#echo "memory used in Gb is "$mem

python randomizeFastq.py paired1.fastq paired2.fastq P
head -n $((4*$numreads)) P1.fastq > paired1.fastq
head -n $((4*$numreads)) P2.fastq > paired2.fastq

echo
echo "Finished subsampling fastq files"
ls

#echo -e 'Starting Assembly: Number of reads is: \t'$(( $( wc -l < $paired1) / 4 ))
#cat $single1 $single2 > single.fastq

# Assembly -m stands for memory allowed
/usr/bin/time -v /local10G/brianyu/tools/Spades-3.9.0/bin/spades.py --only-assembler -k 21,33,55,77 --careful --threads 12 --sc -m 120 -1 paired2.fastq -2 paired1.fastq -o $output_dir

#mv spade_results/contigs.fasta $contigs
#mv spade_results/scaffolds.fasta $scaffolds
#mv spade_results/spades.log $logfile

mv paired* $output_dir

date
echo "Assembly Completed"


# Cleanup 
echo cleanup
date
exit




