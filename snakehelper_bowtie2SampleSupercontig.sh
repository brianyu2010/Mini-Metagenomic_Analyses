#!/bin/bash

###################################
# Using Bowtie2 to align unclustered
# reads from one subsample to the 
# final super contigs obtain through
# a variety of methods. 
# The output is a vector containing
# the coverage depth normalized by 
# contig length.
###################################

scratch=$1
read1=$2
read2=$3
contigs=$4
report=$5
bamfile=$6
bcffile=$7
vcffile=$8
tool_dir=$9
code_dir=${10}
length_thresh=${11}   # this is an integer number 
thread=${12}

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
cp $code_dir/threshold_scaffolds.py $scratch

# filter out scaffolds shorter than contig_thresh      
python threshold_scaffolds.py $length_thresh $contigs temp.fasta

# Build Bowtie2 index
samfile="alignResults"
$tool_dir"bowtie2-2.2.6/bowtie2-build" -f temp.fasta spadeContigs

# Align back with Bowtie2
$tool_dir"bowtie2-2.2.6/bowtie2" --very-sensitive-local -I 0 -X 1000 -p $thread -t -x spadeContigs \
-1 $read1 -2 $read2 -S $samfile.sam 

# Do pileup for each contig
# 2016.08.18 Changed to use bowtie2-2.2.6 and samtools 1.3. I think samtools sort uses -o option now
$tool_dir"samtools-1.3/samtools" view -b -o $samfile".bam" $samfile".sam"
$tool_dir"samtools-1.3/samtools" sort -o $samfile"_sorted.bam" $samfile".bam"
$tool_dir"samtools-1.3/samtools" index $samfile"_sorted.bam"
# in future samtool 1.3 version, -D is replaced with -t DB
$tool_dir"samtools-1.3/samtools" mpileup -f temp.fasta -o $samfile".pile" $samfile"_sorted.bam"
$tool_dir"samtools-1.3/samtools" mpileup -g -t DP -f temp.fasta -o $samfile".bcf" $samfile"_sorted.bam"
$tool_dir"bcftools-1.3/bcftools" call -c -v -o $samfile".vcf" $samfile".bcf"

# Tabulate contig coverage normalized by depth (could use python)
# this statement allows me to sum all the numbers in column 4 corresponding to each id
# in column 1. The output is contig id \t total depth of coverage (not normalized by contig length)
#
# There is no need to sort the first column
# sort -k1 $samfile".pile" > $samfile"_sorted.pile"

# so now use $samfile".pile" directly
awk 'a1==$1 {a2+=$4; next} {print a1,"\t", a2; a1=$1; a2=$4} END {print a1,"\t",a2}' $samfile".pile" > $report

mv $samfile"_sorted.bam" $bamfile
mv $samfile".bcf" $bcffile
mv $samfile".vcf" $vcffile

echo "Alignment Completed"

date
exit


