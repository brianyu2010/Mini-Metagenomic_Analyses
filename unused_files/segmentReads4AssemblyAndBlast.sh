#!/bin/bash

# File Description:
# Calls kmerTrim with various parameters to do assembly processing
# This function can also be used to test paramter space.
# This script needs to be called on scigc-submit-1
#
# Revision History:
# 2014.12.12 Created
# 2014.12.28 Adapted for /local10G/ single cell cluster
# 2015.04.13 Adapted to add Blast
#

# Import variables
#AssemblySGECode=$1	# include full path

##################################
# User Set Directories
##################################

EXP_ROOT=$1
SEEDFILE=$2					# include full path
CODE_DIR=$EXP_ROOT"code_analysis"
RESULT_DIR=$EXP_ROOT"results"
REPORT_DIR=$EXP_ROOT"alignment_reports"

# Check seedfile
if [ ! -e $SEEDFILE ]
then
  echo Seed file $SEEDFILE file  
  exit
fi


	for fastqSize in 100000 #80000 200000 #400000 600000 800000 1000000 1500000 40000
	do

	for kmerSize in 5
	do

	for similarity in 0.98
	do

		echo -e $line'\t'$similarity'\t'$kmerSize'\t'$fastqSize
		# for debuging
		# sbatch -o $REPORT_DIR"kmer_Trim_%A_%a.out" --array=1-27 $CODE_DIR/kmerTrim_spade_v1.sh $SEEDFILE $similarity $kmerSize $fastqSize 
		sbatch -o $REPORT_DIR/"kmerTrimBlast_%A_%a.out" --array=1-6 $CODE_DIR/kmerTrim_spade_blast_v2.sh $SEEDFILE $similarity $kmerSize $fastqSize 

	done

	done

	done


