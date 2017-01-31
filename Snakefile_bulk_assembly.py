###############################################
# Snakemake rules associated with assembly
# blast, quast of bulk subsamples.
# this file must be included into another 
# Snakemake file
#
# 2017.01.13 Added metaspade assembly in addition
#            megahit assembly
###############################################

# rules

rule megahit_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2015.11.30 This rule is still using splitNcluster outputs.
  #            These need to be changed to the P1 P2 S1 S2 from trim
  # 2015.12.01 Changes were made to use trimmed reads directly
  # 2017.01.13 Updated with megahit functions for Bulk
  input: 
    "{subsample}/P1.{subsample}.fastq",
    "{subsample}/P2.{subsample}.fastq",
    "{subsample}/S1.{subsample}.fastq", 
    "{subsample}/S2.{subsample}.fastq"
  output: 
    "{subsample}/megahit_contigs.{subsample}.fasta" 
  params: 
    name="bulk_megahit_assembly", 
    partition=parameters.ix['subsample_assembly_partition','entry'], 
    mem=parameters.ix['subsample_assembly_memory','entry']
  threads: int(parameters.ix['subsample_assembly_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    # perform file size check
    file_size_check = 1;
    for f in input:
      if os.path.getsize(f) == 0:
        file_size_check = 0
    # if file size check passes then proceed
    if file_size_check:
      cp_to_scratch(input,scratch)
      # Assembly using megahit (for metagenomes). 
      # 2015.12.01 Edited to use trimmed reads directly
      # but first delete the output directory so that megahit doesn't complain
      # shell("if [ -d {work_directory}/{wildcards.subsample}/megahit_output_{wildcards.subsample} ]; then rm -rf {work_directory}/{wildcards.subsample}/megahit_output_{wildcards.subsample}; fi") 
      # now do assembly, megahit_output_folder is kept on scratch
      # shell("bash {code_dir}/snakehelper_megahitAssembly.sh {scratch} {input_on_scratch} {resources_dir} {params.mem} {threads} {work_directory}/{wildcards.subsample}/megahit_output_{wildcards.subsample}")
      shell("bash {code_dir}/snakehelper_megahitAssembly.sh {scratch} {input_on_scratch} {resources_dir} {params.mem} {threads} megahit_output_{wildcards.subsample}")
      shell("pwd") # This is to check where I am
      # megahit is the prefix to contig file so OUT_DIR/megahit.contigs.fa
      # shell("cp {wildcards.subsample}/megahit_output_{wildcards.subsample}/megahit.contigs.fa {output[0]}")
      # Now just copy the megahit.contigs.fa file back from scratch megahit_output folder
      shell("cp {scratch}/megahit_output_{wildcards.subsample}/megahit.contigs.fa {output[0]}")
      assert(file_empty(output)),"The contig file is empty."
    else:
      print('Input files have size 0')


rule metaSPAdes_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2017.01.13 Added this to compare with megahit assembly on bulk
  input: 
    "{subsample}/P1.{subsample}.fastq",
    "{subsample}/P2.{subsample}.fastq"
    # "{subsample}/S1.{subsample}.fastq", 
    # "{subsample}/S2.{subsample}.fastq"
  output: 
    "{subsample}/metaSPAdes_contigs.{subsample}.fasta" 
  params: 
    name="bulk_metaSPAdes_assembly", 
    partition=parameters.ix['biosample_assembly_partition','entry'], 
    mem=parameters.ix['biosample_assembly_memory','entry']
  threads: int(parameters.ix['biosample_assembly_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    # perform file size check
    file_size_check = 1;
    for f in input:
      if os.path.getsize(f) == 0:
        file_size_check = 0
    # if file size check passes then proceed
    if file_size_check:
      cp_to_scratch(input,scratch)
      # Assembly using metaSPAdes (for shotgun metagenomes). 
      # now do assembly, keep output folder in scratch
      shell("bash {code_dir}/snakehelper_metaSpadesAssembly.sh {scratch} {input_on_scratch} {code_dir} {tool_dir} {params.mem} {threads} metaSPAdes_output_{wildcards.subsample}")
      shell("pwd") # This is to check where I am
      # megahit is the prefix to contig file so OUT_DIR/megahit.contigs.fa
      shell("cp {scratch}/metaSPAdes_output_{wildcards.subsample}/contigs.fasta {output[0]}")
      assert(file_empty(output)),"The contig file is empty."
    else:
      print('Input files have size 0')


rule align_to_assembly:
  # The output of this rule can be made more informative by including coverage of each contigs.
  input: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq", 
    "{subsample}/contigs.{subsample}.fasta"
  output: "{subsample}/contigCoverage.{subsample}.cnt"
  params: 
    name="Bowtie2Align", 
    partition=parameters.ix['subsample_bowtie2_partition','entry'], 
    mem=parameters.ix['subsample_bowtie2_memory','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Bowtie2 Alignment of reads back to contigs
    shell("bash {code_dir}/snakehelper_bowtie2align2contig.sh {scratch} {input_on_scratch} \
      {output_on_scratch} {tool_dir} {threads} {wildcards.subsample}")
    cp_from_scratch(output, scratch)


rule metaSPAdes_quast:
  input: 
    "{folder}/metaSPAdes_contigs.{k}.fasta"
  output: 
    "{folder}/metaSPAdes_quast_report.{k}.txt"
  params: 
    name="metaSPAdes_quast", 
    partition="general", 
    mem="5000"
  threads: 1
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    output_on_scratch = names_on_scratch(output, scratch)
    input_on_scratch = names_on_scratch(input, scratch)
    cp_to_scratch(input, scratch)
    # Performing Quast
    shell("bash {code_dir}/snakehelper_quast.sh {scratch} {input_on_scratch} {output_on_scratch} {tool_dir} {threads}")
    cp_from_scratch(output, scratch)


rule megahit_quast:
  input: 
    "{folder}/megahit_contigs.{k}.fasta"
  output: 
    "{folder}/megahit_quast_report.{k}.txt"
  params: 
    name="megahit_quast", 
    partition="general", 
    mem="5000"
  threads: 1
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    output_on_scratch = names_on_scratch(output, scratch)
    input_on_scratch = names_on_scratch(input, scratch)
    cp_to_scratch(input, scratch)
    # Performing Quast
    shell("bash {code_dir}/snakehelper_quast.sh {scratch} {input_on_scratch} {output_on_scratch} {tool_dir} {threads}")
    cp_from_scratch(output, scratch)


rule subsample_BLAST:
  input: "{subsample}/contigs.{subsample}.fasta"
  output: "{subsample}/BlastResults.{subsample}.txt"
  params: 
    name="subsample_blast_results", 
    partition=parameters.ix['blast_partition','entry'], 
    mem=parameters.ix['blast_memory','entry'],
    contig_thresh=parameters.ix['subsample_contig_thresh','entry']
  threads: int(parameters.ix['blast_threads','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Blast contigs
    shell("bash {code_dir}/snakehelper_blast.sh {scratch} {input_on_scratch} {output_on_scratch} {tool_dir} {code_dir} {threads} {params.contig_thresh}")
    assert(os.path.isfile(output_on_scratch[0])),"Blast results file does not exist."
    # if file_empty() returns 0 it actually means file is empty
    if not file_empty([output_on_scratch[0]]):
        print('Blast results is empty.')
    #assert(file_empty([output_on_scratch[0]])),"Blast results is empty."
    cp_from_scratch(output, scratch)


