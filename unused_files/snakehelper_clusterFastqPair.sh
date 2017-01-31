#!/bin/bash

#################################
# Fastq Preprocessing 
# Clustering Reads
#################################

scratch=$1
fastq1=$2
fastq2=$3
out1=$4
out2=$5
code_dir=$6
tool_dir=$7
similarity=$8
kmer_len=$9

# if more than 10 arguments use ${10}

echo $scratch
cd $scratch
ls
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
scp $code_dir/ProcessClusteredFastq.py $scratch


# This is for debugging purposes only by reducing file size
# mv $fastq1 temp1.fastq
# mv $fastq2 temp2.fastq
# head -n 1000 temp1.fastq > $fastq1
# head -n 1000 temp2.fastq > $fastq2

echo "Converting First Fastq File"
echo $(( $( wc -l < $fastq1 ) / 4 ))
$tool_dir"fastx/fastq_to_fasta" -Q33 -i $fastq1 -o R1.fasta
head R1.fasta
grep '>' R1.fasta | wc -l

echo "Clustering First Fastq File" $fasta1
/usr/bin/time -v $tool_dir"dnaclust/dnaclust" R1.fasta -l -s $similarity -k $kmer_len > R1clust.txt
# $tool_dir"dnaclust/dnaclust" R1.fasta -l -s $similarity -k $kmer_len > R1clust.txt
date

echo "Converting Second Fastq File"
echo $(( $( wc -l < $fastq2 ) / 4 ))
$tool_dir"fastx/fastq_to_fasta" -Q33 -i $fastq2 -o R2.fasta
head R2.fasta
grep '>' R2.fasta | wc -l

echo "Clustering Second Fastq File" $fasta2
/usr/bin/time -v $tool_dir"dnaclust/dnaclust" R2.fasta -l -s $similarity -k $kmer_len > R2clust.txt
# $tool_dir"dnaclust/dnaclust" R2.fasta -l -s $similarity -k $kmer_len > R2clust.txt
date

# wc -l R*clust.txt
python ProcessClusteredFastq.py R1clust.txt R2clust.txt $fastq1 $fastq2 $out1 $out2
date

echo "Clustering 2 fastq files completed"
echo $(( $( wc -l < $out1 ) / 4 ))
echo $(( $( wc -l < $out2 ) / 4 ))

exit

