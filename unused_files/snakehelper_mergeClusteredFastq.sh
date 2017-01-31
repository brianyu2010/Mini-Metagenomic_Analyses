#!/bin/bash

#################################
# Merging clustered fastq files
# clusterR1_x*.fastq and
# clusterR2_x*.fastq 
# are hard coded. Input file names 
# must be the same.
#################################

scratch=$1
fastq1=$2
fastq2=$3

echo $scratch $fastq1 $fastq2
cd $scratch
ls
date
echo 

# The input filenames are hard coded!!!
ls clusterR1_x*.fastq | sort > R1names.txt
head R1names.txt
echo

ls clusterR2_x*.fastq | sort > R2names.txt 
head R2names.txt
echo

touch $fastq1 $fastq2
while read line; do cat $line >> $fastq1; done < R1names.txt 
while read line; do cat $line >> $fastq2; done < R2names.txt

echo "Fastq merging completed"
exit

