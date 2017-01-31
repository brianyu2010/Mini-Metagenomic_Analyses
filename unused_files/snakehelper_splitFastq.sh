#!/bin/bash

#################################
# Dividing the fastq file into
# small chuncks. This creates files
# with names xaaaab xaaaab etc.
#################################

scratch=$1
CODE_DIR=$2
input1=$3
input2=$4
splitSize=$5

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
scp $CODE_DIR/randomizeFastq.py $scratch

# Sort fastqs randomly
python randomizeFastq.py $input1 $input2 R

# split fastq1 it into differet files
split -l $splitSize R1.fastq -a 5
fastqList=$( ls xaa* | sed -e 's/ /\n/g' )
ls xaa* | xargs -I filename mv filename "R1_"filename".fastq"

echo
# split it into differet files
split -l $splitSize R2.fastq -a 5
ls xaa* | xargs -I filename mv filename "R2_"filename".fastq"
echo 

# Check all fastq files have the same length
for fastqName in ${fastqList[@]}
do
	#echo $fastqName
	if [ $( wc -l < "R1_"$fastqName".fastq" ) -eq $( wc -l < "R2_"$fastqName".fastq" ) ]
	then
		echo 'Fastq files have the same length ' "R1_"$fastqName".fastq" "R2_"$fastqName".fastq"
	else
		echo 'Fastq files have different length ' "R1_"$fastqName".fastq" "R2_"$fastqName".fastq"
	fi
done

echo "Splitting shell script completed"
date

exit

