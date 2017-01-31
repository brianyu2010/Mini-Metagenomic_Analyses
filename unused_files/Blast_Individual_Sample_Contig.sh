#!/bin/bash

#SBATCH --job-name=Blast_Contig_Indiv 	# job name
#SBATCH --ntasks=1  								# number of nodes to allocate per job
#SBATCH --cpus-per-task=12					# cpus (threads) per task
#SBATCH --partition=unrestricted		# partition (queue) to use.
#SBATCH --mem-per-cpu=5333					# memory request per cpu (threads) in MB
##SBATCH --time=72:00:00							# wall clock time (d-H:MM:SS or minutes)
##SBATCH --output=/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/alignment_reports/Trim_Blast_Output.txt
		# this must be a file name

## This is a comment and not recognized by the cluster

## script assumes length of seedfile provided by command line with "-t 1-<file length>"


# Revision History
# 2014.12.18 Created
# 2014.12.25 Updated for Trim and Blast
# 2015.01.24 Updated to Blast assemblies contigs for each sample

date
hostname

##################################
# User Set Directories
##################################

DATA_DIR=/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/141202_M00361_0233_000000000-AB717/
CODE_DIR=/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/code_analysis/


# Tools Directories
tool_dir=/local10G/brianyu/tools/
REF_DIR=/datastore/brianyu/genome_index/
export BLASTDB="/local10G/brianyu/tools/ncbi-blast-2.2.30/db/"
source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/

SEEDFILE=$1	# tab delimited


# Setup Process (User Set Variables)
echo $LOCAL_SATA
if [ -f $SEEDFILE ] # if it's a seedfile
then
	temp=$(head -n $SLURM_ARRAY_TASK_ID $SEEDFILE | tail -n 1)
else	# It's just a folder
	temp=$SEEDFILE
fi
echo $temp
# The cut delimiter is for MiSeq only
sample_id=$(echo $temp | cut -d- -f 2,3)
exp_dir=$DATA_DIR$temp
echo $exp_dir

# Creating Local Copy of Files
echo "create local file copy"
cd $exp_dir/spades_0.98_5_100000
scp scaffolds.fasta $LOCAL_SATA
# scp *.fastq.gz $LOCAL_SATA

cd $LOCAL_SATA

date
echo

##########
# Blast
##########

if [ -s scaffolds.fasta ]
then
	# std is 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
	echo "Starting Blastn"
	echo "Number of Reads is "$(( $( wc -l < scaffolds.fasta ) / 2 ))
	## /usr/bin/time -v 
	$tool_dir"ncbi-blast-2.2.30/bin/blastn" -db "nt" -query scaffolds.fasta -out BlastResults.txt -outfmt '6 std staxids sscinames scomnames sskingdoms stitle qlen'
else
	# Cleanup 
	echo cleanup
	date
	ls
	exit
fi


###########################
# Process Blast Output
###########################

# 15th field is common names, 16th field is kingdoms
if [ -e BlastResults.txt ]
then
	cut -f 15 BlastResults.txt | sort | uniq -c > blastCommonNames.txt
	cut -f 16 BlastResults.txt | sort | uniq -c > blastKindoms.txt
else
	# Cleanup 
	echo cleanup
	date
	ls
	exit
fi

# Copy Back Results
mkdir blast_persample_contigs
mv *.txt blast_persample_contigs
scp -r blast_persample_contigs/ $exp_dir

# Cleanup 
echo cleanup
date
exit




