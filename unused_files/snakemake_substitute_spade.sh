#!/bin/bash

#SBATCH --job-name=Assemble_sub 	# job name
#SBATCH --ntasks=1  								# number of nodes to allocate per job
#SBATCH --cpus-per-task=4						# cpus (threads) per task
#SBATCH --partition=long		# partition (queue) to use.
##SBATCH --mem-per-cpu=5333					# memory request per cpu (threads) in MB
#SBATCH --mem=20000									# total memory
##SBATCH --time=72:00:00							# wall clock time (d-H:MM:SS or minutes)

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
# 2015.05.22 Made into a hack to rescew snakemake assembly section

date
hostname


##################################
# User Set Directories
##################################

# first make a file called assembly_seedfile with folder names and then 
# pass in root directory on local10G
root_dir=$1
folders_file=$1"/assembly_seedfile"

# Tools Directories
tool_dir=/local10G/brianyu/tools/
source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/

num_cores=4
tot_mem=20

# copy files over
cd $LOCAL_SATA
foldername=$( head -n $SLURM_ARRAY_TASK_ID $folders_file | tail -n 1 )
cp $root_dir/$foldername/P*.*.fastq $LOCAL_SATA
cp $root_dir/$foldername/S*.*.fastq $LOCAL_SATA


##################################
# Assembly
##################################

echo -e 'Starting Assembly: Number of reads is: \t'$(( $( wc -l < P1.$foldername.fastq) / 4 ))

cat S1.$foldername.fastq S2.$foldername.fastq > single.fastq

# Assembly -m stands for memory allowed
/usr/bin/time -v python $tool_dir"Spades-3.0.0/bin/spades.py" -k 33,55,77,99,127 --careful --sc -t $num_cores -m $tot_mem \
-1 P1.$foldername.fastq -2 P2.$foldername.fastq -s single.fastq -o spade_output_$foldername > spadeLog.txt
date
echo "Assembly Completed"


# Copy back results
echo "copy back results"

cp spade_output_$foldername/contigs.fasta $root_dir/$foldername/contigs.$foldername.fasta
cp spade_output_$foldername/scaffolds.fasta $root_dir/$foldername/scaffolds.$foldername.fasta
cp -r spade_output_$foldername/ $root_dir/$foldername/


# Cleanup 
echo cleanup
date
exit




