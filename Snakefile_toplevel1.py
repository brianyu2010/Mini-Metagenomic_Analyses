# This snakemake file performs the following actions on one biosample.
# This biosample may contain multiple subsamples and be sequenced multiple times.
# Further, each sequencing run can be missing some subsamples.
# This snakemake file does not perform optimization of thresholding perameters.
#
# 1. quality trimming
# 2. binning and removing overrepresented reads
# 3. assembly
# 4. assembly assessment and realignment
# 5. blast to nt
# 6. combine from different subsamples and assess total results
# 
# This is how to call snakemake
#
# module load use.singlecell
# module load python/3.4.2
# snakemake -n -s snakefile argument (for running just the names)
# snakemake --cluster "sbatch --job-name={params.name} --ntasks=1 --cpus-per-task={threads} \
#           --partition={params.partition} --mem={params.mem} -o slurm_output/%j-{params.name}" -p -j -k
# to kill snakemake process call killall snakemake
#
# Cool Commands
# nohup: doesn't kill a process if a shell command is closed
# git add -A: to update changed, added, and deleted file before committing
#
# Revision History
# 2015.05.18 Brian Yu Created
# 2015.05.26 Updated to include combined analysis as well
# 2015.07.13 Some sample have too many reads for biosample_assembly to run on
#            singlecell cluster. They need to be ran on sherlock.stanford.edu
#            The scripts are modified so that it stops at subsample contig 
#            related computations and mergs corrected reads for input.
# 2015.08.31 Mark added a bigmem node with 128 and 256G and 12 nodes on the cluster
#            The toplevel1.py script is updated to include the toplevel2.py script.
#            Now you only have to run toplevel1.py. Another variable added is work_directory

# Import packages
import os, glob, subprocess
import pandas as pd
import numpy as np
from collections import defaultdict
from snakemake.utils import read_job_properties

# Importing variables NEED TO CHANGE THIS ARGUMENT
# "/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/"
root_folder = config["location"] 

# Importing relavent bio/sub-sample folders, IDs and variables
# The parameters file must be called parameter.txt
parameters = pd.read_table(root_folder+"/code_analysis/parameter.txt", index_col=0)

sample_table = pd.read_table(root_folder+parameters.ix["subsample_location_file",'entry'], header=0, index_col=0)
subsampleIDs = list(set(sample_table.index))

# Pulling out variables from parameter.txt
biosample = parameters.ix["biosample_name",'entry']
code_dir = parameters.ix["code_directory",'entry']
tool_dir = parameters.ix['tool_directory','entry']

# similarity = parameters.ix['similarity','entry']
# kmer_len = parameters.ix['kmer_len','entry']
# split_size = parameters.ix['split_size','entry'] # this is in terms of reads
# contig_thresh = parameters.ix['contig_thresh','entry'] # this is contig length threshold for Blast

"""
# Compute the number of files subcontigs from each subsamples (and total) needs to be split into
depth_table = pd.read_table(root_folder+parameters.ix["subsample_count_file",'entry'], index_col=0, header=None)
# total reads should be the 0th column with header 1. there should only be one column
depth_table.rename(columns={1: 'total_reads'}, inplace=True) 
depth_table['subsample_clust_file_numbers'] = [int(max(1,np.floor(x/int(parameters.ix['subsample_split_size','entry'])))) \
  for x in depth_table['total_reads']]
# biosample_clust_file_numbers = 4
# 2016.08.18 This clustering of reads is no longer used.
biosample_clust_file_numbers = int(max(1, depth_table.sum(0)['total_reads']/int(parameters.ix['biosample_split_size','entry'])))
#print(depth_table)
#print(biosample_clust_file_numbers)
"""

# 2016.02.20 Added this line to avoid multi-line shell commands under
# the new snakemake syntax. This removes errors in lines with &&\
# -e handles nonzero exit status, -u handles unset variables --> -o fail
shell.prefix("set -euo pipefail;")

# Add include files or other snakefile rule files
include: "Snakefile.utils_Mark"
include: "Snakefile.utils_Felix"
include: "Snakefile_helper_Brian.py"
include: "Snakefile_import.py"
include: "Snakefile_subsample_assembly.py"
include: "Snakefile_biosample_assembly.py" # Merging corrected reads and biosample assembly
include: "Snakefile_miniMeta_assembly.py"
# include: "Snakefile_combined_analysis.py" # None of the rules in this files are still used
include: "Snakefile_superContigAnalysis.py" # Aligning subsample reads to supercontigs
# include: "Snakefile_subsampleGenomeSize.py" # Align subsample reads to supercontigs, tabulate genome size recovered. Consolidated

# User defined constants
workdir: "/local10G/brianyu/snakemake_results/"+parameters.ix["biosample_name",0]
work_directory = "/local10G/brianyu/snakemake_results/"+parameters.ix["biosample_name",0]
resources_dir = "/local10G/resources/"

#################################
# A list of all the rules
#################################

rule all:
  # sample should be a list of the subsample names in a biosample. 
  # These are the long names in Miseq runs but ILxxxx-N7xx-N5xx in Nextseq runs
  # input:  expand("{subsample}/BlastResults.{subsample}.txt", subsample=subsampleIDs)
  input: 
    # These are possible outputs to request
    # expand("Combined_Analysis/P1_corrected.{id}.fastq.gz", id=biosample),
    # expand("Combined_Analysis/P2_corrected.{id}.fastq.gz", id=biosample),
    # expand("Combined_Analysis/S_corrected.{id}.fastq.gz", id=biosample),
    # expand("Combined_Analysis/BlastResults.{id}.txt", id=biosample),
    # expand("Combined_Analysis/quast_report.{id}.txt", id=biosample),
    expand("Combined_Analysis/quast_report_miniMeta.{id}.txt", id=biosample),
    # expand("Combined_Analysis/super_contigs.{id}.fasta", id=biosample),
    expand("Combined_Analysis/super_contigs.{id}.alignment_report.txt", id=biosample),
    # expand("Combined_Analysis/subsample_variants.{id}.vcf", id=biosample),
    # expand("{subsample}/contigCoverage.{subsample}.cnt", subsample=subsampleIDs), 
    # expand("{subsample}/BlastResults.{subsample}.txt", subsample=subsampleIDs),
    expand("{subsample}/quast_report.{subsample}.txt", subsample=subsampleIDs),
    # fastqc of fastq reads in subsamples
    expand("{subsample}/P1.{subsample}.fastqc_results.txt", subsample=subsampleIDs),
    expand("{subsample}/P2.{subsample}.fastqc_results.txt", subsample=subsampleIDs),
    # expand("Combined_Analysis/super_contigs.{id}.subsampleGenomeSize.txt", id=biosample)
  params: 
    name="top_level_assembly", 
    partition="general", 
    mem="3000" 
  threads: 1
  version: "1.0"


