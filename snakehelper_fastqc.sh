#!/bin/bash

############################################
# Checking for Overrepresented Sequences
############################################

scratch=$1
input=$2        # format of the file name is {type}.{subsample}.fastq
output=$3
tool_dir=$4

echo $scratch
cd $scratch
date

# Fastqc

# replaces the 6th last character in the name (which is a .) with a _
# The "c" at the end is to change the folder name to *.fastqc
foldername=$( echo $input | rev | sed s/./_/6 | rev )"c"
echo $foldername

$tool_dir"fastqc/fastqc" -j /usr/java/latest/bin/java $input 
sed -n "/Overrepresented/,/END_MODULE/p" $foldername/fastqc_data.txt | head -n -1 | sed '1d' > fastq_overrepresented.txt
ls
if [ -s fastq1_overrepresented.txt ]
then
  sed '1d' fastq1_overrepresented.txt > fastq1_overrepresented2.txt
  mv fastq1_overrepresented2.txt fastq1_overrepresented.txt
  echo "Trimmed fastq contains overrepresented sequences"
else
  echo "Trimmed fastq does not contain overrepresented sequences"
fi

mv $foldername/fastqc_data.txt $output



