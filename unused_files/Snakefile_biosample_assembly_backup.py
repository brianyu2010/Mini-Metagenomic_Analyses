##################################################################
# Snakemake rules associated with 
# Assembly of the entire run, all sequenced combined
# Snakemake file
#
# Revision History:
# 2015.07.13 Added rule biosample_merge_corrected_fastq.
#            This rule will need the spade corrected reads.
# 2017.01.12 Added rule to create YAML file with all corrected
#            fastq.gz files. This file is used by biosample
#            assembly. untrusted-contigs also included
# 2017.02.01 Rmoved a lot of functions that are no longer used.
#            Change some of the names of output files so that
#            they are more consistent with other rules in
#            minimeta contigs
##################################################################

rule create_biosample_assembly_YAML:
  input:
    expand("{subsample}/P1_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/P2_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/S_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/contigs.{subsample}.fasta", subsample=subsampleIDs)
  output:
    # This is in the root work_directory
    "corrected_reads.yaml"
  params:
    name="create_biosample_assembly_YAML",
    partition="general",
    mem="1000" # don't change this
  threads: 1
  version: "1.0"
  run:
    # creating a text file to pass to SPAdes
    numfiles = int(len(input)/4)
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
      # Write untrusted contigs
      d = f.write('  {\n')
      d = f.write('    type: untrusted-contigs,\n')
      d = f.write('    single reads: [\n')
      d = f.write('      ' + ',\n      '.join(input[(3*numfiles):(4*numfiles)]) + '\n')
      d = f.write('    ]\n')
      d = f.write('  }\n')
      d = f.write(']')

"""
rule biosample_merge_corrected_paired_fastq:
  input:
    expand("{subsample}/P1_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/P2_corrected.{subsample}.fastq.gz", subsample=subsampleIDs)
  output:
    "Combined_Analysis/P1_corrected.{id}.fastq.gz",
    "Combined_Analysis/P2_corrected.{id}.fastq.gz"
  params:
    name="biosample_merge_corrected_paired_fastq",
    partition="long",
    mem="63000" # don't change this
  threads: 3
  version: "1.0"
  run:
    # modify input files names with no .gz ending
    input_nozip = modify_zip_file_names(input)
    output_nozip = modify_zip_file_names(output)
    # unzip all the corrected fastq files
    unzip_files(input)
    # The merge statment checks for input and output formates 
    # merge_fastq_files can also differentiate P2 and P2 files
    merge_fastq_files(input_nozip, output_nozip)
    # zip all the correct fastq files back
    zip_back_files(input_nozip)
    zip_back_files(output_nozip)
    # copy output files to /datastore/brianyu/...
    # shell("rsync -avrP {output[0]} {root_folder}/results/")
    # shell("rsync -avrP {output[1]} {root_folder}/results/")


rule biosample_merge_corrected_single_fastq:
  input:
    expand("{subsample}/S_corrected.{subsample}.fastq.gz", subsample=subsampleIDs)
  output:
    "Combined_Analysis/S_corrected.{id}.fastq.gz"
  params:
    name="biosample_merge_corrected_single_fastq",
    partition="long",
    mem="10000" # don't change this
  threads: 1
  version: "1.0"
  run:
    # modify input files names with no .gz ending
    input_nozip = modify_zip_file_names(input)
    output_nozip = modify_zip_file_names(output)
    # unzip all the corrected fastq files
    unzip_files(input)
    # The merge statment checks for input and output formates 
    # merge_fastq_files can also differentiate P2 and P2 files
    shell("cat {input_nozip} > {output_nozip}")
    # zip all the correct fastq files back
    zip_back_files(input_nozip)
    zip_back_files(output_nozip)
    # copy output files to /datastore/brianyu/...
    # shell("rsync -avrP {output} {root_folder}/results/")


# This rule is no longer used 2015.07.15
rule biosample_cluster_reads:
  input: 
    "Combined_Analysis/biosample_P1_{filler}.fastq", 
    "Combined_Analysis/biosample_P2_{filler}.fastq"
  output:
    temp("Combined_Analysis/biosample_cluster_P1_{filler}.fastq"),
    temp("Combined_Analysis/biosample_cluster_P2_{filler}.fastq")
  params:  
    name="biosample_cluster_fastq", 
    partition=parameters.ix['biosample_clust_partition','entry'],
    mem=parameters.ix['biosample_clust_memory','entry'],
    similarity=parameters.ix['biosample_similarity','entry'],
    kmer_len=parameters.ix['kmer_len','entry'] 
  threads: 1
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Clustering
    assert(check_lines(input[0],input[1])),"Input files have different number of lines."
    shell("bash {code_dir}/snakehelper_clusterFastqPair.sh {scratch} {input_on_scratch} {output_on_scratch} {code_dir} {tool_dir} {params.similarity} {params.kmer_len}")
    assert(file_empty(output_on_scratch)),"One or more of the output files are empty."
    assert(check_lines(output_on_scratch[0],output_on_scratch[1])),"Output files have different number of lines."
    assert(check_fastq_ids(output_on_scratch[0],output_on_scratch[1])),"Output fastq files have different read id orders."
    cp_from_scratch(output, scratch)


# This rule is no longer used 2015.07.15
rule biosample_merge_fastq:
  input:
    expand("Combined_Analysis/biosample_cluster_P1_{split}.fastq", split=range(biosample_clust_file_numbers)),
    expand("Combined_Analysis/biosample_cluster_P2_{split}.fastq", split=range(biosample_clust_file_numbers))
  output:
    temp("Combined_Analysis/P1.{id}.fastq"),
    temp("Combined_Analysis/P2.{id}.fastq")
  params:
    name="biosample_merge_clustered_fastq",
    partition="general",
    mem="2000" # don't change this
  threads: 1
  version: "1.0"
  run:
    # The merge statment checks for input and output formates 
    merge_fastq_files(input, output)
"""


rule biosample_spade_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2017.01.12 input is edited to use YAML file and combined subsample contigs
  input: 
    # "Combined_Analysis/P1_corrected.{id}.fastq.gz",
    # "Combined_Analysis/P2_corrected.{id}.fastq.gz",
    # "Combined_Analysis/S_corrected.{id}.fastq.gz"
    "corrected_reads.yaml" # in the root work_directory because of spades paths
  output: 
    "Combined_Analysis/contigs.{id}.fasta", 
    "Combined_Analysis/scaffolds.{id}.fasta"
  params:    
    name="biosample_spade_assembly", 
    partition=parameters.ix['biosample_assembly_partition','entry'], 
    mem=parameters.ix['biosample_assembly_memory','entry']
  threads: int(parameters.ix['biosample_assembly_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    # input_on_scratch = names_on_scratch(input, scratch)
    # output_on_scratch = names_on_scratch(output, scratch)
    # This cp_to_scratch line should not be needed since I don't use those files
    # cp_to_scratch(input, scratch)
    # So far, this is the only rule that uses the variable work_directory 2015.08.31
    print(work_directory)
    # Assembly using spades V3.5.0 much better kmer normalization
    # One single read fastq.gz file passed in. Output directory is set to the "Combined Analysis" Directory
    # If processes breaks, use the conitue flag
    # shell("bash {code_dir}/snakehelper_bigmemAssembly.sh {scratch} {input} {tool_dir} {params.mem} {threads} {work_directory}/Combined_Analysis/spade_output_{wildcards.id}")
    shell("bash {code_dir}/snakehelper_bigmemYAMLAssembly.sh {scratch} {input} {tool_dir} {params.mem} {threads} {work_directory}/Combined_Analysis/spade_output_{wildcards.id}")
    shell("cp Combined_Analysis/spade_output_{wildcards.id}/contigs.fasta {output[0]}")
    shell("cp Combined_Analysis/spade_output_{wildcards.id}/scaffolds.fasta {output[1]}")
    assert(file_empty(output)),"Either the contig or scaffold file is empty."



# rename each super contig with super contig code
rule make_superContigs:
  input: "{folder}/contigs.{id}.fasta"
  output: "{folder}/super_contigs.{id}.fasta", "{folder}/super_contigs_name.{id}.txt"
  params:
    name="make_super_contigs",
    partition="general",
    mem="1000"
  threads: 1
  version: "1.0"
  run: 
    # The set +u; prevents unbound variable errors 2016.08.18
    shell("set +u; source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/ &&\
      python {code_dir}/snakehelper_construct_superContigs.py {input} -o {output[0]} -l {output[1]}")
    assert(file_empty(output)),"Either the supercontig or supercontig names are empty."



rule biosample_BLAST:
  input: "Combined_Analysis/super_contigs.{id}.fasta"
  output: "Combined_Analysis/BlastResults.{id}.txt"
  params:
    name="biosample_blast_results",
    partition="long",
    mem="63000",
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: 12
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Blast contigs
    assert(file_empty([input[0]])),"Input Contig file is empty."
    shell("bash {code_dir}/snakehelper_blast.sh {scratch} {input_on_scratch} {output_on_scratch} {tool_dir} {code_dir} {threads} {params.contig_thresh}")
    assert(file_empty([output_on_scratch[0]])),"Blast output file is empty."
    cp_from_scratch(output, scratch)



"""
rule sspace_scaffold:
  input: 
    expand("{subsample}/P1.{subsample}.fastq", subsample=subsampleIDs), 
    expand("{subsample}/P2.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/S1.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/S2.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/contigs.{subsample}.fasta", subsample=subsampleIDs)
  output: 
    "Combined_Analysis/contigs_sspace.{id}.fasta", 
    "Combined_Analysis/contigs_sspace_names.{id}.fasta",
    "Combined_Analysis/P1.{id}.fastq",
    "Combined_Analysis/P2.{id}.fastq",
    "Combined_Analysis/S1.{id}.fastq",
    "Combined_Analysis/S2.{id}.fastq",
    "Combined_Analysis/library.{id}.txt"
  params:    
    name="sspace_scaffolding", 
    partition="long", 
    mem="63000"
  threads: 12
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    # Cat all read files and contig files
    p1 = [x for x in input if 'P1' in x]
    p2 = [x for x in input if 'P2' in x]
    s1 = [x for x in input if 'S1' in x]
    s2 = [x for x in input if 'S2' in x]
    contigs = [x for x in input if 'contigs' in x]
    print(contigs)
    shell("cat {p1} > {output_on_scratch[2]}")
    shell("cat {p2} > {output_on_scratch[3]}")
    shell("cat {s1} > {output_on_scratch[4]}")
    shell("cat {s2} > {output_on_scratch[5]}")
    shell("source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/ &&\
      python {code_dir}/snakehelper_combine_subsample_contigs.py {contigs} -o {output_on_scratch[0]} -l {output_on_scratch[1]}")
    # use sspace to try to extend contigs
    shell("echo 'Lib1 bowtie {output_on_scratch[2]} {output_on_scratch[3]} 300 0.3 FR' > {output_on_scratch[6]}")
    #shell("echo 'Lib2 bowtie {output_on_scratch[4]}' >> {output_on_scratch[6]}")
    #shell("echo 'Lib3 bowtie {output_on_scratch[5]}' >> {output_on_scratch[6]}")
    shell("head {output_on_scratch[6]}")
    # you don't need the absolute path for this prefix when using sspace
    shell("/local10G/brianyu/tools/sspace-3.0/SSPACE_Standard_v3.0.pl -l {output_on_scratch[6]} \
      -s {output_on_scratch[0]} -k5 -a0.5 -n100 -z{contig_thresh} -x1 -m50 -o10 -r0.8 -b sspace_{wildcards.id}")
    shell("cp -r {scratch}/sspace_{wildcards.id} Combined_Analysis/")
    cp_from_scratch(output, scratch)
"""
