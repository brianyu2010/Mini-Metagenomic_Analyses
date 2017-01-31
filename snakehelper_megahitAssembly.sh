#!/bin/bash

############################################
# Megahit assembly
# Does not necessarily run on the bigmem node
#
# 2017.01.13 Updated names to be megahit instead
#            of spades_xxx 
############################################

# All input should be zipped files

scratch=$1
paired1=$2
paired2=$3
single1=$4
single2=$5
tool_dir=$6    # this is actually /local10G/resources/ for megahit
tot_mem=$7     # this is in MB, need to add 6 zeros
threads=$8
megahit_output_dir=$9

echo $scratch
cd $scratch
date
ls
echo

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
mem=$( echo $tot_mem )000000
echo "memory used in Mb is "$tot_mem

# Here i'm just assuming that the paired end fastq files passed in are not .gz files
echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( cat $paired1 | wc -l ) / 4 ))

if [ -e $single1 ] && [ -e $single2 ]
then
  # zcat handles zip files directly, cat cannot                                      
  cat $single1 $single2 > single.fastq
  echo -e 'Starting Assembly: Number of single end reads is: \t'$(( $( wc -l < single.fastq ) / 4 ))
else
  echo "No single end fastq reads detected"
fi

# Checking if both single files are given
if [ -e $single1 ] && ! [ -e $single2 ]
then
  echo "Only "$single1" found"
fi

if [ -e $single2 ] && ! [ -e $single1 ]
then
  echo "Only "$single2" found"
fi

# Assembly -m stands for memory allowed
# Check if there are single end read files.
# 
# The problem with megahit is that it cannot overwrite the previous OUT_DIR Therefore, you must
# delete the output directory if you want to rerun it. If you want to rerun with different memory
# you still have to delete the folder. Therefore, there's really no need to keep the continue option.
# What I will do is delete the folder in the rule that calls this script.
#
# could also use --k-max to specify the max kmer <=127 That way megahit uses more kmers
# 
# 2017.01.13 Before was using --k-list 33,55,77,99, now use --presets meta (or meta-sensitive)

if [ -e single.fastq ]
then
  if [ -d $megahit_output_dir ]
  then
    /usr/bin/time -v $tool_dir"megahit/megahit" --k-list 21,31,41,51,61,71 -t $threads -m $mem --mem-flag 1 --continue -1 $paired1 -2 $paired2 -r single.fastq -o $megahit_output_dir --out-prefix megahit
  else
    /usr/bin/time -v $tool_dir"megahit/megahit" --k-list 21,31,41,51,61,71 -t $threads -m $mem --mem-flag 1 -1 $paired1 -2 $paired2 -r single.fastq -o $megahit_output_dir --out-prefix megahit
  fi

else

  if [ -d $megahit_output_dir ]
  then
    /usr/bin/time -v $tool_dir"megahit/megahit" --k-list 21,31,41,51,61,71 -t $threads -m $mem --mem-flag 1 --continue -1 $paired1 -2 $paired2 -o $megahit_output_dir --out-prefix megahit
  else
    /usr/bin/time -v $tool_dir"megahit/megahit" --k-list 21,31,41,51,61,71 -t $threads -m $mem --mem-flag 1 -1 $paired1 -2 $paired2 -o $megahit_output_dir --out-prefix megahit
  fi
fi

ls | wc -l
date
echo "Assembly Completed"

exit


