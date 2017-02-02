#!/bin/bash

############################################
# Given a folder containing fastq files
# this script copies those fastq files
# over and then concatenates them together
# while keeping the paired end reads order
############################################

fastqdir=$1
scratch=$2
read1name=$3
read2name=$4

# import files
cd $fastqdir
ls
# You must do this be cause the full path somehow doesn't work.
cp *.fastq.gz $scratch
cd $scratch

# unzip files
pwd
ls
ls *.fastq.gz > zipped_file_names.txt
while read line
do
  gzip -d $line
done < zipped_file_names.txt

# sort names for cat
ls *_R1_001.fastq | sort > read1files.txt
head read1files.txt
ls *_R2_001.fastq | sort > read2files.txt
head read2files.txt

# cat files and then remove them
# touch $read1name
while read line
do
  cat $line >> $read1name
  echo -e "Number of reads in read1.fastq is: "$(( $( wc -l < $line ) / 4 ))
  rm $line
done < read1files.txt


# touch $read2name
while read line
do
  cat $line >> $read2name
  echo -e "Number of reads in read2.fastq is: "$(( $( wc -l < $line ) / 4 ))
  rm $line
done < read2files.txt

echo
echo -e "Number of reads in Read1.output is: "$(( $( wc -l < $read1name) / 4 ))
echo -e "Number of reads in Read2.output is: "$(( $( wc -l < $read2name) / 4 ))
echo 

echo "Combining fastq files completed"
date
exit


