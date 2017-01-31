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

# Import packages
import os, glob, subprocess
import pandas as pd
from collections import defaultdict
from snakemake.utils import read_job_properties

# Importing variables NEED TO CHANGE THIS ARGUMENT
# "/datastore/brianyu/2014.11.21_YNP_LowerGeyserBasin/"
root_folder = config["location"] 

# Importing relavent bio/sub-sample folders, IDs and variables
sample_table = pd.read_table(root_folder+"/code_analysis/subsamples.txt", header=0, index_col=0)
subsampleIDs = list(set(sample_table.index))
parameters = pd.read_table(root_folder+"/code_analysis/parameter.txt", index_col=0)
assert("biosample_name" in parameters.index)
biosample = parameters.ix["biosample_name",0]
assert("code_directory" in parameters.index)
code_dir = parameters.ix["code_directory",0]
assert("tool_directory" in parameters.index)
tool_dir = parameters.ix['tool_directory',0]
assert("similarity" in parameters.index)
similarity = parameters.ix['similarity',0]
assert("kmer_len" in parameters.index)
kmer_len = parameters.ix['kmer_len',0]
assert("split_size" in parameters.index)
split_size = parameters.ix['split_size',0]
assert('contig_thresh' in parameters.index)
contig_thresh = parameters.ix['contig_thresh',0] # this is contig length threshold for Blast

# Add include files or other snakefile rule files
include: "Snakefile.utils_Mark"
include: "Snakefile.utils_Felix"
include: "Snakefile.import"
include: "Snakefile.subsample_assembly"
include: "Snakefile.combined_analysis"

# User defined constants
workdir: "/local10G/brianyu/snakemake_results/"+parameters.ix["biosample_name",0]


#################################
# A list of all the rules
#################################

rule combine:
  input: 
    #expand("Combined_Analysis/quast_report.{id}.txt", id=biosample),
    #expand("Combined_Analysis/BlastResults.{id}.txt", id=biosample)
    #expand("Combined_Analysis/super_contigs_distribution.{id}.txt", id=biosample),
    #expand("Combined_Analysis/super_contigs.{id}.similarity_matrix.txt", id=biosample),
    #expand("Combined_Analysis/subsample_contigs.{id}.similarity_matrix.txt", id=biosample),
    #expand("Combined_Analysis/subsample_species_abundance.{id}.txt", id=biosample)
    #expand("Combined_Analysis/contigs_sspace.{id}.fasta", id=biosample)
  params: 
    name="combined_assembly", 
    partition="general", 
    mem="5000" 
  threads: 1
  version: "1.0"

