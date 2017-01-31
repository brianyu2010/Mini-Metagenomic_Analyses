#!/bin/bash

#SBATCH --job-name=Test_Assemble 	# job name
#SBATCH --ntasks=1           # number of nodes to allocate per job
#SBATCH --cpus-per-task=12		# cpus (threads) per task
#SBATCH --partition=long		# partition (queue) to use.
##SBATCH --mem-per-cpu=5333		# memory request per cpu (threads) in MB
#SBATCH --mem=63000			# total memory
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
# 2017.01.12 Only for testing different spade or other assembly methods
# 2017.01.12 Added untrusted and trusted contigs

date
hostname

inputfile=$1
output_dir=$2
contigs=$3

scratch=$LOCAL_SATA

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/

# Assembly -m stands for memory allowed
/usr/bin/time -v /local10G/brianyu/tools/Spades-3.9.0/bin/spades.py --only-assembler -k 21,33,55,77 --threads 12 -m 63 --dataset $inputfile --trusted-contigs $contigs -o $output_dir

date
echo "Assembly Completed"


# Cleanup 
echo cleanup
date
exit




