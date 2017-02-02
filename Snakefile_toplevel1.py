# This snakemake file performs the following actions on one biosample.
# This biosample may contain multiple subsamples and be sequenced multiple times.
# Further, each sequencing run can be missing some subsamples.
# This snakemake file does not perform optimization of thresholding perameters.
#
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
# 2017.01.31 Updated for sherlock.stanford.edu, using hns partition for bigmem jobs

# Import packages
import os, glob, subprocess
import pandas as pd
import numpy as np
from collections import defaultdict
from snakemake.utils import read_job_properties

# Importing variables NEED TO CHANGE THIS ARGUMENT
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
python2_env = parameters.ix['python2_environment_name','entry']
python3_env = parameters.ix['python3_environment_name','entry']
work_directory_base = parameters.ix['work_directory_base','entry']
work_directory = work_directory_base+'/'+parameters.ix["biosample_name",'entry']

# 2016.02.20 Added this line to avoid multi-line shell commands under
# the new snakemake syntax. This removes errors in lines with &&\
# -e handles nonzero exit status, -u handles unset variables --> -o fail
# shell.prefix("set -euo pipefail;")

# Add include files or other snakefile rule files
# include: "Snakefile.utils_Mark"
include: "Snakefile.utils_Felix"
include: "Snakefile_helper_Brian.py"
include: "Snakefile_import.py"
include: "Snakefile_subsample_assembly.py"
include: "Snakefile_biosample_assembly.py" # Merging corrected reads and biosample assembly
include: "Snakefile_miniMeta_assembly.py"
include: "Snakefile_superContigAnalysis.py" # Aligning subsample reads to supercontigs

# User defined constants
workdir: work_directory
# workdir: "/scratch/users/brianyu/"+parameters.ix["biosample_name",0]
# resources_dir = "/local10G/resources/"

#################################
# A list of all the rules
#################################

rule all:
  # sample should be a list of the subsample names in a biosample. 
  # These are the long names in Miseq runs but ILxxxx-N7xx-N5xx in Nextseq runs
  # input:  expand("{subsample}/BlastResults.{subsample}.txt", subsample=subsampleIDs)
  input: 
    # These are possible outputs to request
    expand("Combined_Analysis/super_contigs.{id}.fasta", id=biosample),
    expand("Combined_Analysis/super_contigs.{id}.alignment_report.txt", id=biosample),
    expand("Combined_Analysis/subsample_variants.{id}.vcf", id=biosample),
    expand("{subsample}/quast_report.{subsample}.txt", subsample=subsampleIDs),
    expand("{subsample}/P1.{subsample}.fastqc_results.txt", subsample=subsampleIDs),
    expand("{subsample}/P2.{subsample}.fastqc_results.txt", subsample=subsampleIDs),
    expand("Combined_Analysis/super_contigs.{id}.subsampleGenomeSize.txt", id=biosample)
  params: 
    name="top_level_assembly",
    qos="normal",
    time="30:00",
    partition="normal", 
    mem="4000" 
  threads: 1
  version: "2.0"


