#!/bin/bash

###################################
# Blast
###################################

scratch=$1
contigs=$2
report=$3
tool_dir=$4
code_dir=$5
num_cores=$6       # this is an integer number
length_thresh=$7   # this is an integer number

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
export BLASTDB="/local10G/brianyu/tools/ncbi-blast-2.2.30/db/"
cp $code_dir/threshold_scaffolds.py $scratch

# filter out scaffolds shorter than contig_thresh
python threshold_scaffolds.py $length_thresh $contigs temp.fasta

##############################################
# std is 'qseqid sseqid pident length mismatch 
# gapopen qstart qend sstart send evalue bitscore'
# You should compare bitscore 
##############################################

echo "Starting Blastn"
echo "Number of Reads is "$(( $( wc -l < temp.fasta ) / 2 ))

## /usr/bin/time -v 
# returns 1 alignment and uses 8 threads
echo "number of cores used is "$num_cores
$tool_dir"ncbi-blast-2.2.30/bin/blastn" -db "nt" -query temp.fasta -num_threads $num_cores -num_alignments 1 -out $report -outfmt '6 qseqid sseqid bitscore pident length evalue staxids sscinames scomnames sskingdoms stitle qlen'
echo "Blast Completed"

date
exit


