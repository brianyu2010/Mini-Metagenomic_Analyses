#!/bin/bash

#SBATCH --job-name=Trim_And_Blast 	# job name
#SBATCH --ntasks=1  								# number of nodes to allocate per job
#SBATCH --cpus-per-task=12					# cpus (threads) per task
#SBATCH --partition=unrestricted		# partition (queue) to use.
#SBATCH --mem-per-cpu=5333					# memory request per cpu (threads) in MB
##SBATCH --time=72:00:00							# wall clock time (d-H:MM:SS or minutes)
##SBATCH --output=/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/alignment_reports/Trim_Blast_Output.txt
		# this must be a file name

## This is a comment and not recognized by the cluster

## script assumes length of seedfile provided by command line with "-t 1-<file length>"


# Revision History
# 2014.12.18 Created
# 2014.12.25 Updated for Trim and Blast

date
hostname

##################################
# User Set Directories
##################################

DATA_DIR=/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/141202_M00361_0233_000000000-AB717/
CODE_DIR=/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/code_analysis/


# Tools Directories
tool_dir=/local10G/brianyu/tools/
REF_DIR=/datastore/brianyu/genome_index/
export BLASTDB="/local10G/brianyu/tools/ncbi-blast-2.2.30/db/"
source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/

SEEDFILE=$1	# tab delimited


# Setup Process (User Set Variables)
echo $LOCAL_SATA
if [ -f $SEEDFILE ] # if it's a seedfile
then
	temp=$(head -n $SLURM_ARRAY_TASK_ID $SEEDFILE | tail -n 1)
else	# It's just a folder
	temp=$SEEDFILE
fi
echo $temp
# The cut delimiter is for MiSeq only
sample_id=$(echo $temp | cut -d- -f 2,3)
exp_dir=$DATA_DIR$temp
echo $exp_dir

# Creating Local Copy of Files
echo "create local file copy"
cd $exp_dir
scp *.fastq.gz $LOCAL_SATA
scp $tool_dir"trimmomatic/adapters/Combined_PE_V2.fa" $LOCAL_SATA/adapterSeqs.fa

cd $LOCAL_SATA

# Combine Read Files
cat *_R1_001.fastq.gz > R1.fastq.gz
cat *_R2_001.fastq.gz > R2.fastq.gz
gzip -d R1.fastq.gz
gzip -d R2.fastq.gz

date
echo

# Trimmomatic variables
default_seedMismatch=3
default_palen_th=30
default_minAdapterLen=3
default_slidingWindow=10
default_slidingQual=25
default_maxInfoLen=120
default_maxInfo_th=0.3

default_option_string="ILLUMINACLIP:adapterSeqs.fa:$default_seedMismatch:$default_palen_th:10:$default_minAdapterLen:TRUE SLIDINGWINDOW:$default_slidingWindow:$default_slidingQual MAXINFO:$default_maxInfoLen:$default_maxInfo_th LEADING:30 TRAILING:30 MINLEN:30"

#####################################################
# Use SeqPrep to join paired end reads when possible
#####################################################

$tool_dir"seqprep/SeqPrep" -f R1.fastq -r R2.fastq -1 R1single.fastq.gz -2 R2single.fastq.gz -s merged.fastq.gz
gzip -d *.fastq.gz

###################
# Trimming Reads
###################

if [ -s R1single.fastq ] && [ -s R2single.fastq ] && [ -s merged.fastq ]
then
	echo "Start Trimming with Trimmomatic"
	# /usr/java/latest/bin/java
	/usr/java/latest/bin/java -jar $tool_dir"trimmomatic/trimmomatic-0.30.jar" PE -phred33 -trimlog pairtrim.log R1single.fastq R2single.fastq pairedout1.fastq singleout1.fastq pairedout2.fastq singleout2.fastq $default_option_string 
	/usr/java/latest/bin/java -jar $tool_dir"trimmomatic/trimmomatic-0.30.jar" SE -phred33 -trimlog mergetrim.log merged.fastq mergedout.fastq $default_option_string 
	echo "Trimming Completed"
else
	# Cleanup 
	echo cleanup
	date
	ls
	exit
fi

############################################
# Checking for Overrepresented Sequences
############################################

if [ -s pairedout1.fastq ]
then
$tool_dir"fastqc/fastqc" -j /usr/java/latest/bin/java pairedout1.fastq
sed -n "/Overrepresented/,/END_MODULE/p" pairedout1_fastqc/fastqc_data.txt | head -n -1 | sed '1d' > fastq1_overrepresented.txt
if [ -s fastq1_overrepresented.txt ]
then
  sed '1d' fastq1_overrepresented.txt > fastq1_overrepresented2.txt
  mv fastq1_overrepresented2.txt fastq1_overrepresented.txt
  echo "Trimmed forward fastq contains overrepresented sequences"
else
  echo "Trimmed forward fastq does not contain overrepresented sequences"
fi
mv pairedout1_fastqc/fastqc_data.txt pairedout1_fastqc.txt
fi

if [ -s pairedout2.fastq ]
then
$tool_dir"fastqc/fastqc" -j /usr/java/latest/bin/java pairedout2.fastq
sed -n "/Overrepresented/,/END_MODULE/p" pairedout2_fastqc/fastqc_data.txt | head -n -1 | sed '1d' > fastq2_overrepresented.txt
if [ -s fastq2_overrepresented.txt ]
then
  sed '1d' fastq2_overrepresented.txt > fastq2_overrepresented2.txt
  mv fastq2_overrepresented2.txt fastq2_overrepresented.txt
  echo "Trimmed reverse fastq contains overrepresented sequences"
else
  echo "Trimmed reverse fastq does not contain overrepresented sequences"
fi
mv pairedout2_fastqc/fastqc_data.txt pairedout2_fastqc.txt
fi

if [ -s mergedout.fastq ]
then
$tool_dir"fastqc/fastqc" -j /usr/java/latest/bin/java mergedout.fastq
sed -n "/Overrepresented/,/END_MODULE/p" mergedout_fastqc/fastqc_data.txt | head -n -1 | sed '1d' > mergedout_overrepresented.txt
if [ -s mergedout_overrepresented.txt ]
then
  sed '1d' mergedout_overrepresented.txt > mergedout_overrepresented2.txt
  mv mergedout_overrepresented2.txt mergedout_overrepresented.txt
  echo "Trimmed merged fastq contains overrepresented sequences"
else
  echo "Trimmed merged fastq does not contain overrepresented sequences"
fi
mv mergedout_fastqc/fastqc_data.txt mergedout_fastqc.txt
fi

#############################
# Convert Fastq to Fasta
#############################

if [ -s pairedout1.fastq ] && [ -s pairedout2.fastq ] && [ -s mergedout.fastq ] 
	then
	echo "Concatenate All fastq files together"
	cat pairedout1.fastq pairedout2.fastq mergedout.fastq > combinedReads.fastq
	echo "Converting Fastq to Fasta"
	$tool_dir"fastx/fastq_to_fasta" -Q33 -i combinedReads.fastq -o combinedReads.fasta
	date
else
	# Cleanup 
	echo cleanup
	date
	ls
	exit
fi

##########
# Blast
##########

# std is 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
echo "Starting Blastn"
echo "Number of Reads is "$(( $( wc -l < combinedReads.fasta ) / 2 ))
## /usr/bin/time -v 
$tool_dir"ncbi-blast-2.2.30/bin/blastn" -db "nt" -query combinedReads.fasta -out BlastResults.txt -outfmt '6 std staxids sscinames scomnames sskingdoms stitle qlen'

###########################
# Process Blast Output
###########################

# 15th field is common names, 16th field is kingdoms
if [ -e BlastResults.txt ]
then
	cut -f 15 BlastResults.txt | sort | uniq -c > blastCommonNames.txt
	cut -f 16 BlastResults.txt | sort | uniq -c > blastKindoms.txt
else
	# Cleanup 
	echo cleanup
	date
	ls
	exit
fi

# Copy Back Results
mkdir blast_output
mv *.txt blast_output
mv pairedout1.fastq blast_output
mv pairedout2.fastq blast_output
mv mergedout.fastq blast_output
scp -r blast_output/ $exp_dir

# Cleanup 
echo cleanup
date
exit




