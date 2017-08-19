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
# 2015.12.03 Change the structure of this toplevel file to handle bulk metagenomic sequencing
#            samples. What this means is that this toplevel file no longer has the combined
#            analysis portions. In addition, since Spades is not build to do assembly on
#            complex metagenomes, this file uses Megahit for assembly.
#            This file uses rules in Snakefile_import.py and a new file Snakefile_bulk_assembly.py
# 2017.01.31 Updated for sherlock.stanford.edu, using hns partition for bigmem jobs
# 2017.02.01 Added bulk sample processing to proceed at the same time
# 2017.07.08 Added functionality to perform genome reassembly from mini-metagenomic and bulk sub-samples
#            Need to add another entry into parameters file called genome_cluster_file
# 2017.07.21 Made bowtie2 into 1 rule, added checkm and prokka rules

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
parameters = pd.read_table(root_folder+"/code_analysis/parameter-sherlock.txt", index_col=0)

# Getting Minimetagenomic subsamples
sample_table = pd.read_table(root_folder+parameters.ix["subsample_location_file",'entry'], header=0, index_col=0)
subsampleIDs = list(set(sample_table.index))

# Getting bulk sample processing. Bulk samples should be in the same folder as subsamples
bulk_flag = parameters.ix['bulk_samples_present','entry']
if bulk_flag == 'Yes' or bulk_flag == 'yes' or bulk_flag == 'Y' or bulk_flag == 'y':
  bulk_table = pd.read_table(root_folder+parameters.ix['bulksample_location_file','entry'], header=0, index_col=0)
  bulksampleIDs = list(set(bulk_table.index))

# Getting genomes from first round clustering.
genome_table = pd.read_table(root_folder+parameters.ix["genome_cluster_file",'entry'], header=0, index_col=0)
genomeIDs = list(set(genome_table.index))
genomeIDs = ['{:03.0f}'.format(x) for x in genomeIDs]
# print(genomeIDs)
# genomeIDs = range(1,int(parameters.ix["genome_number",'entry']))
# special python 3 str.format function

# Pulling out variables from parameter.txt
biosample = parameters.ix["biosample_name",'entry']
code_dir = parameters.ix["code_directory",'entry']
tool_dir = parameters.ix['tool_directory','entry']
python2_env = parameters.ix['python2_environment_name','entry']
python3_env = parameters.ix['python3_environment_name','entry']
prokka_env = parameters.ix['prokka_environment_name','entry']
work_directory_base = parameters.ix['work_directory_base','entry']
work_directory = work_directory_base+'/'+parameters.ix["biosample_name",'entry']

# Add include files or other snakefile rule files
# include: "Snakefile.utils_Mark"
include: "Snakefile.utils_Felix"
include: "Snakefile_helper_Brian.py"

# User defined constants
workdir: work_directory

#################################
# A list of all the rules
#################################

if bulk_flag =='Yes' or bulk_flag == 'yes' or bulk_flag == 'Y'or bulk_flag =='y':
  rule all_withbulk:
    input:
      expand("Genome_Reassembly/genome_contigs_miniMeta.{genome}.fasta", genome=genomeIDs),
      expand("Genome_Reassembly/genome_contigs_withBulk.{genome}.fasta", genome=genomeIDs),
      expand("Genome_Reassembly/prokka_annotation_{genome}/genome_annotation.{genome}.txt", genome=genomeIDs),
      "Genome_Reassembly/checkm_lineage_wf_completed.txt"
    params:
      name="top_level_genome_reassembly",
      qos="normal",
      time="10:00",
      partition="quake,normal,owner",
      mem="4000"
    threads: 1
    version: "2.0"
else:
  rule all_withbulk:
    input:
      expand("Genome_Reassembly/genome_contigs_miniMeta.{genome}.fasta", genome=genomeIDs),
      # expand("Genome_Reassembly/genome_contigs_withBulk.{genome}.fasta", genome=genomeIDs)
    params:
      name="top_level_genome_reassembly",
      qos="normal",
      time="10:00",
      partition="quake,normal,owner",
      mem="4000"
    threads: 1
    version: "2.0"


rule make_genome_index:
  input:
    "Genome_groups/genome_cluster.{genome}.fasta"
  output:
    temp("Genome_groups/genome_cluster_{genome}.1.bt2l"),
    temp("Genome_groups/genome_cluster_{genome}.2.bt2l"),
    temp("Genome_groups/genome_cluster_{genome}.3.bt2l"),
    temp("Genome_groups/genome_cluster_{genome}.4.bt2l"),
    temp("Genome_groups/genome_cluster_{genome}.rev.1.bt2l"),
    temp("Genome_groups/genome_cluster_{genome}.rev.2.bt2l")
  params:
    name="make_genome_index",
    qos="normal",
    time="12:00:00",
    partition="quake,normal",
    mem="64000"
  threads: 5
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing contig trimming and bowtie2-build
    # The variable work_directory is used in this part
    shell("""
      basename=$(echo {output[0]} | cut -d/ -f2 | cut -d. -f1)
      echo $basename
      date; cd {scratch}; 
      source activate {python3_env}
      bowtie2-build --quiet --large-index --threads {threads} -f {input_on_scratch} $basename
      source deactivate; date
      """)
    # Move the rest of the reads back
    cp_from_scratch(output, scratch)


rule extract_reads_bulk:
  input:
    "{bulksample}/P1.{bulksample}.fastq",
    "{bulksample}/P2.{bulksample}.fastq",
    "{bulksample}/S1.{bulksample}.fastq",
    "{bulksample}/S2.{bulksample}.fastq",
    "Genome_groups/genome_cluster_{genome}.1.bt2l",
    "Genome_groups/genome_cluster_{genome}.2.bt2l",
    "Genome_groups/genome_cluster_{genome}.3.bt2l",
    "Genome_groups/genome_cluster_{genome}.4.bt2l",
    "Genome_groups/genome_cluster_{genome}.rev.1.bt2l",
    "Genome_groups/genome_cluster_{genome}.rev.2.bt2l"
  output: # once works, make output files temp
    temp("Genome_groups/bulkReads_{bulksample}_GenomeCluster_{genome}_P1.fastq"),
    temp("Genome_groups/bulkReads_{bulksample}_GenomeCluster_{genome}_P2.fastq"),
    temp("Genome_groups/bulkReads_{bulksample}_GenomeCluster_{genome}_S.fastq")
  params:
    name="extract_reads_bulk",
    qos="normal",
    time=parameters.ix['reassembly_bowtie2_time','entry'],
    mem_per_core="6G",
    partition=parameters.ix['reassembly_bowtie2_partition','entry'],
    mem=parameters.ix['reassembly_bowtie2_memory','entry']
  threads: int(parameters.ix['reassembly_bowtie2_thread','entry'])
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    # bowtie2 can handle fastq.gz files
    shell("""
      date; source activate {python3_env}
      # basename=$(echo {input[4]} | cut -d/ -f2 | cut -d. -f1)
      basename=$(echo {input[4]} | cut -d. -f1)
      echo $basename
      # cd {scratch}
      cat {input[2]} {input[3]} > {scratch}/single_reads.fastq
      echo 'Using local alignment to include as many reads as possible reads'
      bowtie2 --time --phred33 --very-sensitive-local -I 100 -X 2000 -p {threads} --al {scratch}/alignedSingle --al-conc {scratch}/alignedPaired -x $basename -1 {input[0]} -2 {input[1]} -U {scratch}/single_reads.fastq -S {scratch}/alignResults.sam
      mv {scratch}/alignedPaired.1 {output[0]}
      mv {scratch}/alignedPaired.2 {output[1]}
      mv {scratch}/alignedSingle {output[2]}
      source deactivate
      """)


rule extract_reads_minimeta_allSubSamples:
  input:
    expand("{subsample}/P1.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/P2.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/S1.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/S2.{subsample}.fastq", subsample=subsampleIDs),
    "Genome_groups/genome_cluster_{genome}.1.bt2l",
    "Genome_groups/genome_cluster_{genome}.2.bt2l",
    "Genome_groups/genome_cluster_{genome}.3.bt2l",
    "Genome_groups/genome_cluster_{genome}.4.bt2l",
    "Genome_groups/genome_cluster_{genome}.rev.1.bt2l",
    "Genome_groups/genome_cluster_{genome}.rev.2.bt2l"
  output: # output files could be temp
    temp("Genome_groups/miniMetaReads_GenomeCluster_{genome}_P1.fastq"),
    temp("Genome_groups/miniMetaReads_GenomeCluster_{genome}_P2.fastq"),
    temp("Genome_groups/miniMetaReads_GenomeCluster_{genome}_S.fastq")
  params:
    name="extract_reads_minimeta_allSubSamples",
    qos="normal",
    time=parameters.ix['reassembly_bowtie2_time','entry'],
    mem_per_core="6G",
    partition=parameters.ix['reassembly_bowtie2_partition','entry'],
    mem=parameters.ix['reassembly_bowtie2_memory','entry']
  threads: int(parameters.ix['reassembly_bowtie2_thread','entry'])
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    # Performing bowtie2 alignment --time --phred33 --un <path> --un-conc <path>
    numSubSamples = len(subsampleIDs)
    p1 = input[(0 * numSubSamples) : (1 * numSubSamples)]
    p2 = input[(1 * numSubSamples) : (2 * numSubSamples)]
    s1 = input[(2 * numSubSamples) : (3 * numSubSamples)]
    s2 = input[(3 * numSubSamples) : (4 * numSubSamples)]
    single_reads_name = ['' for i in range(numSubSamples)]
    for i in range(numSubSamples):
      single_reads_name[i] = scratch + '/single_reads.' + subsampleIDs[i] + '.fastq'
      shell("cat {s1[i]} {s2[i]} > {single_reads_name[i]}")
    p1 = ','.join(p1)
    p2 = ','.join(p2)
    s = ','.join(single_reads_name)
    shell("""
      date; source activate {python3_env}
      basename=$(echo {input[-1]} | cut -d. -f1)
      echo $basename
      echo 'Using local alignment to remove reads'
      bowtie2 --time --phred33 --very-sensitive-local -I 100 -X 2000 -p {threads} --al {scratch}/alignedSingle --al-conc {scratch}/alignedPaired -x $basename -1 {p1} -2 {p2} -U {s} -S {scratch}/alignResults.sam
      mv {scratch}/alignedPaired.1 {output[0]}
      mv {scratch}/alignedPaired.2 {output[1]}
      mv {scratch}/alignedSingle {output[2]}
      source deactivate
      """)


rule assemble_genome_with_bulk:
  input:
    # expand("Genome_groups/miniMetaReads_{subsample}_GenomeCluster_{{genome}}_P1.fastq", subsample=subsampleIDs),
    # expand("Genome_groups/miniMetaReads_{subsample}_GenomeCluster_{{genome}}_P2.fastq", subsample=subsampleIDs),
    # expand("Genome_groups/miniMetaReads_{subsample}_GenomeCluster_{{genome}}_S.fastq", subsample=subsampleIDs),
    "Genome_groups/miniMetaReads_GenomeCluster_{genome}_P1.fastq",
    "Genome_groups/miniMetaReads_GenomeCluster_{genome}_P2.fastq",
    "Genome_groups/miniMetaReads_GenomeCluster_{genome}_S.fastq",
    expand("Genome_groups/bulkReads_{bulksample}_GenomeCluster_{{genome}}_P1.fastq", bulksample=bulksampleIDs),
    expand("Genome_groups/bulkReads_{bulksample}_GenomeCluster_{{genome}}_P2.fastq", bulksample=bulksampleIDs),
    expand("Genome_groups/bulkReads_{bulksample}_GenomeCluster_{{genome}}_S.fastq", bulksample=bulksampleIDs),
    "Genome_groups/genome_cluster.{genome}.fasta"
  output:
    "genome_reassembly_withBulk.{genome}.yaml",
    "Genome_Reassembly/genome_contigs_withBulk.{genome}.fasta",
    "Genome_Reassembly/genome_scaffolds_withBulk.{genome}.fasta",
    "Genome_Reassembly/quast_genome_withBulk.{genome}.txt"
  params:
    name="assemble_genome_with_bulk",
    qos="normal",
    time=parameters.ix['reassembly_spades_time','entry'],
    partition=parameters.ix['reassembly_spades_partition','entry'],
    mem=parameters.ix['reassembly_spades_memory','entry'], #"128000",
    contig_thresh=parameters.ix['reassembly_contig_thresh','entry'],
    kmer=parameters.ix['biosample_spades_kmerlist','entry'] # "21,33,44,77,99"
  threads: int(parameters.ix['reassembly_spades_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    # creating a text file to pass to SPAdes
    # numfiles_subsample = len(subsampleIDs)
    numfiles_subsample = 1
    numfiles_bulksample = len(bulksampleIDs)
    print(numfiles_subsample)
    print(numfiles_bulksample)
    with open(output[0],'w') as f:
      # Write file names of left reads
      d = f.write('[\n')
      d = f.write('  {\n')
      d = f.write('    orientation: fr,\n')
      d = f.write('    type: paired-end,\n')
      d = f.write('    left reads: [\n')
      d = f.write('      ' + ',\n      '.join(input[0:numfiles_subsample]) + ',\n')
      d = f.write('      ' + ',\n      '.join(input[(3*numfiles_subsample):(3*numfiles_subsample+numfiles_bulksample)]) + '\n')
      d = f.write('    ],\n')
      # Write file names of right reads
      d = f.write('    right reads: [\n')
      d = f.write('      ' + ',\n      '.join(input[numfiles_subsample:(2*numfiles_subsample)]) + ',\n')
      d = f.write('      ' + ',\n      '.join(input[(3*numfiles_subsample+numfiles_bulksample):(3*numfiles_subsample+2*numfiles_bulksample)]) + '\n')
      d = f.write('    ]\n')
      d = f.write('  },\n')
      # Write file names of single reads
      d = f.write('  {\n')
      d = f.write('    type: single,\n')
      d = f.write('    single reads: [\n')
      d = f.write('      ' + ',\n      '.join(input[(2*numfiles_subsample):(3*numfiles_subsample)]) + ',\n')
      d = f.write('      ' + ',\n      '.join(input[(3*numfiles_subsample+2*numfiles_bulksample):(3*numfiles_subsample+3*numfiles_bulksample)]) + '\n')
      d = f.write('    ]\n')
      d = f.write('  },\n')
      # Write untrusted contigs for this assembly, there should only be one file
      d = f.write('  {\n')
      d = f.write('    type: untrusted-contigs,\n')
      d = f.write('    single reads: [\n')
      d = f.write('      ' + input[-1] + '\n')
      d = f.write('    ]\n')
      d = f.write('  }\n')
      d = f.write(']')
    print('Completed creating YAML file')
    # Performing spade assembly, not using --meta
    # quast.py --plots-format svg --gene-finding -meta -m {params.contig_thresh} -t {threads} -o {scratch}/quast_output {output[1]}
    shell("""
      echo {scratch}; date; pwd; echo
      source activate {python2_env}
      mem=$( echo {params.mem} | rev | cut -c 4- | rev )
      echo 'memory used in Gb is '$mem
      spades_output_dir={scratch}/spades_Re-Assembly_{wildcards.genome}
      source activate {python3_env}
      time spades.py --phred-offset 33 --sc --careful -k {params.kmer} -t {threads} -m $mem --dataset {output[0]} -o $spades_output_dir
      cp $spades_output_dir/contigs.fasta {output[1]}
      cp $spades_output_dir/scaffolds.fasta {output[2]}
      source activate {python2_env}
      # --meta expired after 400 days. --glimmer is not for metagenomics
      quast.py -m {params.contig_thresh} -t {threads} -o $spades_output_dir/quast_output {output[1]}
      cp $spades_output_dir/quast_output/report.txt {output[3]}
      date; source deactivate
      """)
    assert(file_empty(output)),"Either the contig or scaffold file is empty."


rule assemble_genome_minimeta_only:
  input:
    "Genome_groups/miniMetaReads_GenomeCluster_{genome}_P1.fastq",
    "Genome_groups/miniMetaReads_GenomeCluster_{genome}_P2.fastq",
    "Genome_groups/miniMetaReads_GenomeCluster_{genome}_S.fastq",
    # expand("Genome_groups/miniMetaReads_{subsample}_GenomeCluster_{{genome}}_P1.fastq", subsample=subsampleIDs),
    # expand("Genome_groups/miniMetaReads_{subsample}_GenomeCluster_{{genome}}_P2.fastq", subsample=subsampleIDs),
    # expand("Genome_groups/miniMetaReads_{subsample}_GenomeCluster_{{genome}}_S.fastq", subsample=subsampleIDs),
    "Genome_groups/genome_cluster.{genome}.fasta"
  output:
    "genome_reassembly_miniMeta.{genome}.yaml",
    "Genome_Reassembly/genome_contigs_miniMeta.{genome}.fasta",
    "Genome_Reassembly/genome_scaffolds_miniMeta.{genome}.fasta",
    "Genome_Reassembly/quast_genome_miniMeta.{genome}.txt"
  params:
    name="assemble_genome_miniMeta",
    qos="normal",
    time=parameters.ix['reassembly_spades_time','entry'],
    partition=parameters.ix['reassembly_spades_partition','entry'],
    mem=parameters.ix['reassembly_spades_memory','entry'],
    contig_thresh=parameters.ix['reassembly_contig_thresh','entry'],
    kmer=parameters.ix['biosample_spades_kmerlist','entry'] # "21,33,44,77,99"
  threads: int(parameters.ix['reassembly_spades_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    # creating a text file to pass to SPAdes
    numfiles = int((len(input)-1)/3)
    print(numfiles)
    with open(output[0],'w') as f:
      # Write file names of left reads
      d = f.write('[\n')
      d = f.write('  {\n')
      d = f.write('    orientation: fr,\n')
      d = f.write('    type: paired-end,\n')
      d = f.write('    left reads: [\n')
      d = f.write('      ' + ',\n      '.join(input[0:numfiles]) + '\n')
      d = f.write('    ],\n')
      # Write file names of right reads
      d = f.write('    right reads: [\n')
      d = f.write('      ' + ',\n      '.join(input[numfiles:(2*numfiles)]) + '\n')
      d = f.write('    ]\n')
      d = f.write('  },\n')
      # Write file names of single reads
      d = f.write('  {\n')
      d = f.write('    type: single,\n')
      d = f.write('    single reads: [\n')
      d = f.write('      ' + ',\n      '.join(input[(2*numfiles):(3*numfiles)]) + '\n')
      d = f.write('    ]\n')
      d = f.write('  },\n')
      # Write untrusted contigs for this assembly, there should only be one file
      d = f.write('  {\n')
      d = f.write('    type: untrusted-contigs,\n')
      d = f.write('    single reads: [\n')
      d = f.write('      ' + input[-1] + '\n')
      d = f.write('    ]\n')
      d = f.write('  }\n')
      d = f.write(']')
    print('Completed creating YAML file')
    # Performing spade assembly, not using --meta
    # quast.py --plots-format svg --gene-finding -meta -m {params.contig_thresh} -t {threads} -o {scratch}/quast_output {output[1]}
    shell("""
      echo {scratch}; date; pwd; echo
      source activate {python2_env}
      mem=$( echo {params.mem} | rev | cut -c 4- | rev )
      echo 'memory used in Gb is '$mem
      spades_output_dir={scratch}/spades_Re-Assembly_{wildcards.genome}
      source activate {python3_env}
      time spades.py --phred-offset 33 --sc --careful -k {params.kmer} -t {threads} -m $mem --dataset {output[0]} -o $spades_output_dir
      cp $spades_output_dir/contigs.fasta {output[1]}
      cp $spades_output_dir/scaffolds.fasta {output[2]}
      source activate {python2_env}
      # do not set min contig thresh for this, default is 500bp, --gene-finding --glimmer
      quast.py -t {threads} -o $spades_output_dir/quast_output {output[1]}
      cp $spades_output_dir/quast_output/report.txt {output[3]}
      date; source deactivate
      """)
    assert(file_empty(output)),"Either the contig or scaffold file is empty."



rule prokka_annotate_genome:
  input:
    "Genome_Reassembly/genome_contigs_withBulk.{genome}.fasta"
  output:
    "Genome_Reassembly/prokka_annotation_{genome}/genome_annotation.{genome}.txt"
  params:
    name="prokka_annotate_genome",
    qos="normal",
    time="12:00:00",
    partition="quake,normal,owners",
    mem="64000"
  threads: 5
  version: "1.0"
  run:
    scratch = os.environ["LOCAL_SCRATCH"]
    shell("""
      date; echo
      outdir=$( echo {output} | cut -d/ -f 1-2 )
      prefix=$( echo {output} | cut -d/ -f3 | cut -d. -f 1-2 )
      echo $outdir; echo $prefix; echo
      source activate {prokka_env}
      prokka --metagenome --force --outdir $outdir --prefix $prefix --cpus {threads} {input}
      source deactivate
      """)


rule checkm_genome:
  input:
    # expand("Genome_Reassembly/genome_contigs_withBulk.{genome}.fasta", genome=genomeIDs)
    expand("Genome_Reassembly/genome_scaffolds_withBulk.{genome}.fasta", genome=genomeIDs)
  output:
    "Genome_Reassembly/checkm_lineage_wf_completed.txt"
  params:
    name="checkm_genome",
    qos="normal",
    time="12:00:00",
    partition="quake,normal",
    contig_thresh=parameters.ix['reassembly_contig_thresh','entry'],
    mem="128000"
  threads: 10
  version: "1.0"
  run:
    # Manage Files
    scratch = os.environ["LOCAL_SCRATCH"]
    # input_on_scratch = names_on_scratch(input, scratch)
    # cp_to_scratch(input, scratch)
    # Perform Process
    genome_folder = "genome_bin"
    checkm_result = "checkm_result"
    shell("""
      mkdir {scratch}/{genome_folder}
      date; echo
      # cd {scratch}; ls
      source activate {python2_env}
      for contig_file in {input}
      do
        echo $contig_file
        scratch_contig_file="{scratch}/{genome_folder}/"$( echo $contig_file | rev | cut -d/ -f1 | rev )
        echo $scratch_contig_file
        python {code_dir}/process_scaffolds.py --lengthThresh {params.contig_thresh} $contig_file $scratch_contig_file
      done
      ls {scratch}; echo; ls {scratch}/{genome_folder}; echo; echo
      # mv *.fasta {genome_folder}
      """)
    shell("""
      curdir=$(pwd)
      output_dir=$( echo {output} | cut -d/ -f1 )
      cd {scratch}; date; echo; ls; echo;
      source activate {python2_env}
      mkdir temp_dir
      checkm lineage_wf -t {threads} --pplacer_threads {threads} --extension fasta --tmpdir temp_dir {genome_folder} {checkm_result}
      source deactivate; echo;
      # right now you are still in {scratch}
      echo "Move checkm results back"; echo; ls; echo
      if [ -d $curdir/$output_dir/{checkm_result} ]
      then
        rm -rf $curdir/$output_dir/{checkm_result}
      fi
      mv {checkm_result} $curdir/$output_dir/
      """)
    shell("""
      output_dir=$( echo {output} | cut -d/ -f1 ) 
      if [ -d $output_dir/{checkm_result} ]
      then
        date; echo "checkm Completed"; echo
        touch {output}
      fi
      """)

