#!/bin/bash

#################################
# Quality Trimming Fastq Files
#################################

scratch=$1
in1=$2
in2=$3
paired1=$4
single1=$5
paired2=$6
single2=$7
tool_dir=$8
thread=$9

echo $scratch
cd $scratch
date


# This is for debugging purposes only by reducing file size
# echo "Subset input read files for debugging."
# mv $in1 temp1.fastq
# mv $in2 temp2.fastq
# head -n 10000 temp1.fastq > $in1
# head -n 10000 temp2.fastq > $in2



cp $tool_dir"trimmomatic/adapters/Combined_PE_V2.fa" $scratch/adapterSeqs.fa

# Trimmomatic variables
default_seedMismatch=3
default_palen_th=30
default_minAdapterLen=3
default_slidingWindow=10
default_slidingQual=25
default_maxInfoLen=120
default_maxInfo_th=0.3

default_option_string="ILLUMINACLIP:adapterSeqs.fa:$default_seedMismatch:$default_palen_th:10:$default_minAdapterLen:TRUE SLIDINGWINDOW:$default_slidingWindow:$default_slidingQual MAXINFO:$default_maxInfoLen:$default_maxInfo_th LEADING:30 TRAILING:30 MINLEN:30"

echo "Start Trimming with Trimmomatic"
# /usr/java/latest/bin/java

# output order is P1.fastq S1.fastq P2.fastq S2.fastq
/usr/java/latest/bin/java -jar $tool_dir"trimmomatic/trimmomatic-0.30.jar" PE -phred33 -threads $thread -trimlog pairtrim.log $in1 $in2 $paired1 $single1 $paired2 $single2 $default_option_string 

echo "Trimming Completed"
date
exit

