#!/bin/bash

# Used to submit jobs to singlecell-submit clusters
# Revision History
# 2014.12.25 Brian Yu

reports_folder=/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/alignment_reports/
array_num=$1
script=$2
seedfile=$3

# sbatch -o $reports_folder"Trim_Blast_%A_%a.out" --array=1-$array_num $script $seedfile 

sbatch -o $reports_folder"Blast_Indiv_Contigs_%A_%a.out" --array=1-$array_num $script $seedfile 




