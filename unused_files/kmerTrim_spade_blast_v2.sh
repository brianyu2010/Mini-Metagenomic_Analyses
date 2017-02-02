#!/bin/bash

#SBATCH --job-name=Cluster_n_Assemble 	# job name
#SBATCH --ntasks=1  								# number of nodes to allocate per job
#SBATCH --cpus-per-task=12						# cpus (threads) per task
#SBATCH --partition=long		# partition (queue) to use.
##SBATCH --mem-per-cpu=5333					# memory request per cpu (threads) in MB
#SBATCH --mem=60000									# total memory
##SBATCH --time=72:00:00							# wall clock time (d-H:MM:SS or minutes)

## This is a comment and not recognized by the cluster

## script assumes length of seedfile provided by command line with "--array=1-<file length>"

# File Description:
# This file is run on /local10G/brianyu/ and is called by segmentReads4Assembly.sh

# Revision History
# 2014.05.16 Created
# 2014.10.27 Edited to use DNA clust first and then spade and do quast. Also keep kmer frequency
# 2014.12.11 Added ProcessClusteredFastq.py to re-order paired end reads
# 2014.12.12 Removed seedfile as an entry. Instead, experiment folder, similarity, kmer, and output_folder are passed in
#						 This file now must be called by another shell script which passes in these values
# 2014.12.14 Actually Added back the functionality of passing in a seed file
# 2014.12.28 Adapted for single cell cluster at /local10G/
# 2015.01.07 Changed to Version 2 where no reference is available and included trimming step.
#            Still does not include blast step after assembly
# 2015.04.13 Added blast function and removed Quast. This is for YNP environmental samples.

date
hostname


##################################
# User Set Directories
##################################

EXP_ROOT=/datastore/brianyu/2015.03.19_YNP_gDNACtrl
DATA_DIR=$EXP_ROOT/150406_M00361_0252_000000000-AD92F/
DATA2_DIR=$EXP_ROOT/150402_M00361_0251_000000000-AD8VY/
CODE_DIR=$EXP_ROOT/code_analysis/
RESULT_DIR=$EXP_ROOT/results/
REPORT_DIR=$EXP_ROOT/alignment_reports/

# Tools Directories
tool_dir=/local10G/brianyu/tools/
REF_DIR=/datastore/brianyu/genome_index/
source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
export BLASTDB="/local10G/brianyu/tools/ncbi-blast-2.2.30/db/"

expfolder=$1			# tab delimited, basiclaly seedfile
similarity=$2			# similarity in dnaclust for testing purposes
kmer_len=$3				# kmer length for clustering
splitSize=$4			# number of lines of each smaller fastq
# outputFolder=$5		# output folder the spade_#_# folder will go into
num_cores=12
tot_mem=60


# Setup Process (User Set Variables)
echo $LOCAL_SATA
if [ -f $expfolder ] # if it's a seedfile
then
	temp=$(head -n $SLURM_ARRAY_TASK_ID $expfolder | tail -n 1)
else	# It's just a folder
	temp=$expfolder
fi
echo $temp
# the delimiter is for nextseq sample names
# Will just return N7xx-N5xx
sample_id=$(echo $temp | cut -d_ -f 6 | cut -d- -f 2,3) # Only for YNP control sample names on Miseq
exp_dir=$DATA_DIR$temp
echo $exp_dir

# Intermediate File Names
samfile="alignResults"
ClustPair="ClustPair"
# all_coord="all.coord"
# mapped_list="mapped.list"
# mapped_coord="mapped.coord"
alignment_file="contigCoverage.cnt"



###################################
# Creating Local Copy of Files
###################################

echo "create local file copy"

echo $exp_dir
cd $exp_dir
cp *.fastq.gz $LOCAL_SATA
cd $LOCAL_SATA
cat *_R1_001.fastq.gz > A1.fastq.gz
cat *_R2_001.fastq.gz > A2.fastq.gz
gzip -d A1.fastq.gz
gzip -d A2.fastq.gz
# rm *_R1_001.fastq.gz
# rm *_R2_001.fastq.gz
echo $DATA2_DIR$temp
cd $DATA2_DIR$temp
cp *.fastq.gz $LOCAL_SATA
cd $LOCAL_SATA
cat *_R1_001.fastq.gz > B1.fastq.gz
cat *_R2_001.fastq.gz > B2.fastq.gz
gzip -d B1.fastq.gz
gzip -d B2.fastq.gz
# rm *_R1_001.fastq.gz
# rm *_R2_001.fastq.gz
# combined unzipped fastq
cat A1.fastq B1.fastq > read1.fastq
cat A2.fastq B2.fastq > read2.fastq

# for debugging only
# head -n 100000 A1.fastq > read1.fastq
# head -n 100000 A2.fastq > read2.fastq

rm A*.fastq
rm B*.fastq


#################################
# Creating Local Copy of Scripts
#################################

echo "create local script copy"
cd $CODE_DIR
scp ProcessClusteredFastq.py $LOCAL_SATA
scp randomizeFastq.py $LOCAL_SATA
scp threshold_scaffolds.py $LOCAL_SATA
scp $tool_dir"trimmomatic/adapters/Combined_PE_V2.fa" $LOCAL_SATA/adapterSeqs.fa

cd $LOCAL_SATA
date
echo


# Check Arguments
echo "Checking Argument Values"
echo -e $kmer_len'\t'$similarity'\t'$splitSize'\t'$outputFolder
spadeFolder="spades_"$similarity"_"$kmer_len"_"$splitSize


#################################
# Quality Trimming Fastq Files
#################################

# Trimmomatic variables
default_seedMismatch=3
default_palen_th=30
default_minAdapterLen=3
default_slidingWindow=10
default_slidingQual=25
default_maxInfoLen=120
default_maxInfo_th=0.3

default_option_string="ILLUMINACLIP:adapterSeqs.fa:$default_seedMismatch:$default_palen_th:10:$default_minAdapterLen:TRUE SLIDINGWINDOW:$default_slidingWindow:$default_slidingQual MAXINFO:$default_maxInfoLen:$default_maxInfo_th LEADING:30 TRAILING:30 MINLEN:30"

if [ -s read1.fastq ] && [ -s read2.fastq ]
then
	echo "Start Trimming with Trimmomatic"
	# /usr/java/latest/bin/java
	/usr/java/latest/bin/java -jar $tool_dir"trimmomatic/trimmomatic-0.30.jar" PE -phred33 -trimlog pairtrim.log read1.fastq read2.fastq P1.fastq S1.fastq P2.fastq S2.fastq $default_option_string 
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

if [ -s P1.fastq ]
then
$tool_dir"fastqc/fastqc" -j /usr/java/latest/bin/java P1.fastq
sed -n "/Overrepresented/,/END_MODULE/p" P1_fastqc/fastqc_data.txt | head -n -1 | sed '1d' > fastq1_overrepresented.txt
if [ -s fastq1_overrepresented.txt ]
then
  sed '1d' fastq1_overrepresented.txt > fastq1_overrepresented2.txt
  mv fastq1_overrepresented2.txt fastq1_overrepresented.txt
  echo "Trimmed forward fastq contains overrepresented sequences"
else
  echo "Trimmed forward fastq does not contain overrepresented sequences"
fi
mv P1_fastqc/fastqc_data.txt P1_fastqc.txt
fi

if [ -s P2.fastq ]
then
$tool_dir"fastqc/fastqc" -j /usr/java/latest/bin/java P2.fastq
sed -n "/Overrepresented/,/END_MODULE/p" P2_fastqc/fastqc_data.txt | head -n -1 | sed '1d' > fastq2_overrepresented.txt
if [ -s fastq2_overrepresented.txt ]
then
  sed '1d' fastq2_overrepresented.txt > fastq2_overrepresented2.txt
  mv fastq2_overrepresented2.txt fastq2_overrepresented.txt
  echo "Trimmed reverse fastq contains overrepresented sequences"
else
  echo "Trimmed reverse fastq does not contain overrepresented sequences"
fi
mv P2_fastqc/fastqc_data.txt P2_fastqc.txt
fi


#################################
# Dividing the fastq file into
# small chuncks. This creates files
# with names xaaaab xaaaab etc.
#################################

if [ -s P1.fastq ] && [ -s P2.fastq ]
then

	# Sort fastqs randomly
	python randomizeFastq.py P1.fastq P2.fastq R

	# split fastq1 it into differet files
	split -l $splitSize R1.fastq -a 5
	fastqList=$( ls xaa* | sed -e 's/ /\n/g' )
	ls xaa* | xargs -I filename mv filename "R1_"filename".fastq"

	# split it into differet files
	split -l $splitSize R2.fastq -a 5
	ls xaa* | xargs -I filename mv filename "R2_"filename".fastq"

	# Check all fastq files have the same length
	for fastqName in ${fastqList[@]}
	do
		# echo $fastqName
		if [ $( wc -l < "R1_"$fastqName".fastq" ) -eq $( wc -l < "R2_"$fastqName".fastq" ) ]
		then
			echo 'Fastq files have the same length ' "R1_"$fastqName".fastq" "R2_"$fastqName".fastq"
		else
			echo 'Fastq files have different length ' "R1_"$fastqName".fastq" "R2_"$fastqName".fastq"
		fi
	done

	# Create the two cumulative fastq files
	touch $ClustPair"1.fastq"
	touch $ClustPair"2.fastq"

else

	# Cleanup 
	echo cleanup
	date
	ls
	exit

fi

#################################
# Fastq Preprocessing 
# For All Segments
#################################

for fastqName in ${fastqList[@]}
do

	cp "R1_"$fastqName".fastq" R1.fastq
	cp "R2_"$fastqName".fastq" R2.fastq

	#################################
	# Fastq Preprocessing 
	#################################

	if [ $similarity != 1.00 ]
	then

		echo "Converting First Fastq File"
		$tool_dir"fastx/fastq_to_fasta" -Q33 -i R1.fastq -o R1.fasta
		echo "Clustering First Fastq File" "R1_"$fastqName".fastq"
		# /usr/bin/time -v $tool_dir"dnaclust/dnaclust" R1.fasta -l -s $similarity -k $kmer_len > R1clust.txt
    $tool_dir"dnaclust/dnaclust" R1.fasta -l -s $similarity -k $kmer_len > R1clust.txt
		date
		echo "Converting Second Fastq File"
		$tool_dir"fastx/fastq_to_fasta" -Q33 -i R2.fastq -o R2.fasta
		echo "Clustering Second Fastq File" "R2_"$fastqName".fastq"
		# /usr/bin/time -v $tool_dir"dnaclust/dnaclust" R2.fasta -l -s $similarity -k $kmer_len > R2clust.txt
    $tool_dir"dnaclust/dnaclust" R2.fasta -l -s $similarity -k $kmer_len > R2clust.txt
		date
		# wc -l R*clust.txt
		python ProcessClusteredFastq.py R1clust.txt R2clust.txt R1.fastq R2.fastq temppaired
		date

	else

		mv R1.fastq temppaired1.fastq
		mv R2.fastq temppaired2.fastq

	fi # if [ $similarity != 1.00 ]


	##################################
	# Appending Processed Fastq Files
	#
	# Note that the same read won't 
	# appear twice in the combined list
	##################################

	cat temppaired1.fastq >> $ClustPair"1.fastq"
	cat temppaired2.fastq >> $ClustPair"2.fastq"

done



##################################
# Assembly
##################################

if [ -s $ClustPair"1.fastq" ] && [ -s $ClustPair"2.fastq" ]

then

	echo -e 'Starting Assembly: Number of reads is: \t'$(( $( wc -l < $ClustPair"1.fastq") / 4 ))

	# Assembly -m stands for memory allowed
	/usr/bin/time -v python $tool_dir"Spades-3.0.0/bin/spades.py" -k 55,77,99,127 --careful --sc -m $tot_mem \
	-1 $ClustPair"1.fastq" -2 $ClustPair"2.fastq" -o $spadeFolder > spadeLog.txt
	date
	echo "Assembly Completed"

else

	echo "unable to find fastq files for assembly"
	ls
	# Cleanup 
	echo cleanup
	date
  exit

fi


###################################
# Quantifying Assembly Statistics
###################################

if [ -f $spadeFolder/scaffolds.fasta ]

then

  # Process scaffold stats with Quast
	python $tool_dir"quast-2.3/quast.py" $spadeFolder/scaffolds.fasta -o quast_output > $spadeFolder/quastLog.txt
	date
	echo "Quast Completed"

  # Align back with Bowtie2
	$tool_dir"bowtie2-2.1.0/bowtie2-build" -f $spadeFolder/scaffolds.fasta spadeContigs > $spadeFolder/bowtie2BuildLog.txt
	$tool_dir"bowtie2-2.1.0/bowtie2" --very-sensitive-local -I 0 -X 1000 -p 8 -t -x spadeContigs \
  -1 $ClustPair"1.fastq" -2 $ClustPair"2.fastq" -S $samfile.sam > $spadeFolder/bowtie2AlignLog.txt
	echo -e "New_Sample\t"$sample_id > $alignment_file
	# Get all coordinates
	echo -e "Total_Reads\t"$( $tool_dir"samtools-0.1.19/samtools" view -Sf 0x001 $samfile".sam" | cut -d: -f 4-7 \
	| cut -f 1 | sort | uniq | wc -l ) >> $alignment_file
	# Uniquely mapped coordinates
	echo -e "Uniquely_Mapped_Reads\t"$( $tool_dir"samtools-0.1.19/samtools" view -Sf 0x003 $samfile".sam" | cut -d: -f 4-7 | cut -f 1-9 \
  | sort | uniq | wc -l ) >> $alignment_file
	# All mapped coordinates (even only one of the paired reads
	echo -e "All_Mapped_Reads\t"$( $tool_dir"samtools-0.1.19/samtools" view -Sf 0x001 -F 0x00c $samfile".sam" | cut -d: -f 4-7 | cut -f 1-9 \
  | sort | uniq | wc -l ) >> $alignment_file
	date
	echo "Alignment Completed"

else

	echo "unable to find scaffolds"
  cp spadeLog.txt $exp_dir/
	ls
	# Cleanup 
	echo cleanup
	date
  exit

fi



###################################
# Blast
###################################

if [ -f $spadeFolder/scaffolds.fasta ]

then

  # filter out scaffolds shorter than 2kb
  cd $LOCAL_SATA
  python threshold_scaffolds.py 900 $spadeFolder/scaffolds.fasta scaffolds.fasta

  # std is 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
	echo "Starting Blastn"
	echo "Number of Reads is "$(( $( wc -l < scaffolds.fasta ) / 2 ))
	## /usr/bin/time -v 
  # returns 1 alignment and uses 8 threads
	$tool_dir"ncbi-blast-2.2.30/bin/blastn" -db "nt" -query scaffolds.fasta -num_threads $num_cores -out BlastResults.txt -outfmt '6 std staxids sscinames scomnames sskingdoms stitle qlen'
  echo "Blast Completed"


	###########################
	# Process Blast Output
	###########################

	# 15th field is common names, 16th field is kingdoms
	if [ -e BlastResults.txt ]
	then

		cut -f 15 BlastResults.txt | sort | uniq -c > blastCommonNames.txt
		cut -f 16 BlastResults.txt | sort | uniq -c > blastKindoms.txt


		# Copy back results
		echo "copy back results"

		# first copy the report files into spades output
		mv $ClustPair"1.fastq" $spadeFolder/
		mv $ClustPair"2.fastq" $spadeFolder/
    mv $alignment_file $spadeFolder/
    mv quast_output/report.txt $spadeFolder/quast_report.txt
    mv P1_fastqc.txt $spadeFolder/
    mv P2_fastqc.txt $spadeFolder/
		mv spadeLog.txt $spadeFolder/
		mv BlastResults.txt $spadeFolder/
		mv blastCommonNames.txt $spadeFolder/
		mv blastKindoms.txt $spadeFolder/
    mv $spadeFolder/scaffolds.fasta $spadeFolder/$sample_id"_scaffolds.fasta"
    mv scaffolds.fasta $spadeFolder/$sample_id"_subset_scaffolds.fasta"
		# scp $alignment_file $spadeFolder/
		rm -rf $spadeFolder/K*/
		rm -rf $spadeFolder/tmp/
		rm -rf $spadeFolder/misc/
		# Copy the folder back to output folder
		scp -r $spadeFolder $exp_dir/
		# Now remove the folder
		rm -rf $spadeFolder/

	else

		echo "unable to find scaffolds"
		ls
		# Cleanup 
		echo cleanup
		date
		exit

	fi


else

	echo "unable to find contigs"
	ls
	# Cleanup 
	echo cleanup
	date
  exit

fi


# Cleanup 
echo cleanup
date
exit




