###############################################
# Snakemake rules associated with assembly
# blast, quast of a single subsample.
# this file must be included into another 
# Snakemake file
###############################################

# rules

rule spade_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2015.11.30 This rule is still using splitNcluster outputs.
  #            These need to be changed to the P1 P2 S1 S2 from trim
  # 2015.12.01 Changes were made to use trimmed reads directly
  input: 
    #"{subsample}/ClustPair1.{subsample}.fastq", 
    #"{subsample}/ClustPair2.{subsample}.fastq"
    "{subsample}/P1.{subsample}.fastq",
    "{subsample}/P2.{subsample}.fastq",
    "{subsample}/S1.{subsample}.fastq", 
    "{subsample}/S2.{subsample}.fastq"
  output: 
    "{subsample}/contigs.{subsample}.fasta", 
    "{subsample}/scaffolds.{subsample}.fasta",
    "{subsample}/P1_corrected.{subsample}.fastq.gz",
    "{subsample}/P2_corrected.{subsample}.fastq.gz",
    "{subsample}/S_corrected.{subsample}.fastq.gz",
    "{subsample}/SPAdes_report.{subsample}.log",
    "{subsample}/SPAdes_params.{subsample}.txt"
  params: 
    name="subsample_spade_assembly", 
    partition=parameters.ix['subsample_assembly_partition','entry'], 
    mem=parameters.ix['subsample_assembly_memory','entry']
  threads: int(parameters.ix['subsample_assembly_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    # perform file size check
    file_size_check = 1;
    for f in input:
      if os.path.getsize(f) == 0:
        file_size_check = 0
    # if file size check passes then proceed
    if file_size_check:
      cp_to_scratch(input, scratch)
      # Assembly using spades V3.5. Only include output_on_scratch[0:2] here 
      # because the snakemake_spadeassembly file is not written to handle the 
      # other two outputs. However you can copy those back later from SCRATCH
      # by the way, [0:2] actually only contains 0th and 1st entry
      # shell("bash {code_dir}/snakehelper_spadeAssembly.sh {scratch} {input_on_scratch} fakeR1 fakeR2 {output_on_scratch[0]} {output_on_scratch[1]} {tool_dir} {params.mem} {threads} spade_output_{wildcards.subsample}")
      # 2015.12.01 Edited to use trimmed reads directly
      shell("bash {code_dir}/snakehelper_spadeAssembly.sh {scratch} {input_on_scratch} {output_on_scratch[0]} {output\
_on_scratch[1]} {tool_dir} {params.mem} {threads} spade_output_{wildcards.subsample}") 
      # 2016.09.23 Could think about commenting the following line out.
      # shell("rsync -avrP {scratch}/spade_output_{wildcards.subsample} {wildcards.subsample}/")
      # edit contig names and copy back file
      # first check if the contig and scaffolds are empty
      assert(file_empty(output_on_scratch[0:2])),"Either contig or scaffold files from SPAdes is empty."
      with open(output_on_scratch[0],'r') as infile, open(output[0],'w') as outfile:
          for line in infile:
              if '>' in line:
                  # Basically all the contigs from the beginning will have
                  # the same naming schemes
                  a = outfile.write('>SubSample_'+wildcards.subsample+'_'+line[1:])
              else:
                  a = outfile.write(line)
      # edit scaffold names and copy back file
      with open(output_on_scratch[1],'r') as infile, open(output[1],'w') as outfile:
          for line in infile:
              if '>' in line:
                  a = outfile.write('>SubSample_'+wildcards.subsample+'_'+line[1:])
              else:
                  a = outfile.write(line)
      # Copy corrected P1 and P2 reads from scratch (unpaired reads are left)
      # These two files are in .fastq.gz format after copying
      # The names of these two files depend on the input files
      shell("""cp {scratch}/correctedreads1.fastq.gz {output[2]}
        cp {scratch}/correctedreads2.fastq.gz {output[3]}
        cp {scratch}/correctedreads_unpaired.fastq.gz {output[4]}
        cp {scratch}/spade_output_{wildcards.subsample}/spades.log {output[5]}
        cp {scratch}/spade_output_{wildcards.subsample}/params.txt {output[6]}""")
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


rule general_quast:
  input: "{folder}/contigs.{k}.fasta"
  output: "{folder}/quast_report.{k}.txt"
  params: 
    name="general_quast", 
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


