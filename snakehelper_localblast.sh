#!/bin/bash

###################################
# Local Blast for Snakemake
###################################

scratch=$1
query=$2
contigs=$3
report=$4
tool_dir=$5
code_dir=$6
num_cores=$7       # this is an integer number
length_thresh=$8   # this is an integer number

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
cp $code_dir/threshold_scaffolds.py $scratch

# filter out scaffolds shorter than contig_thresh
python threshold_scaffolds.py $length_thresh $contigs mydb.fasta
python threshold_scaffolds.py $length_thresh $query query.fasta

# Make local blastn database
$tool_dir"ncbi-blast-2.2.30/bin/makeblastdb" -in mydb.fasta -parse_seqids -dbtype nucl

# std is 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
echo "Starting Blastn"
echo "Number of Contigs is "$(( $( wc -l < query.fasta ) / 2 ))

# doing blast of query against mydb
echo "number of cores used is "$num_cores
/usr/bin/time -v $tool_dir"ncbi-blast-2.2.30/bin/blastn" -db mydb.fasta -query query.fasta -num_threads $num_cores -out tempBlastOut.txt -outfmt '6 qseqid sseqid bitscore pident length'

sort -k 1 tempBlastOut.txt > $report

echo "Blast Completed"

date
exit


