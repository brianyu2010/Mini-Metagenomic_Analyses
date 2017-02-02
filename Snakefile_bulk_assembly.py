###############################################
# Snakemake rules associated with assembly
# blast, quast of bulk subsamples.
# this file must be included into another 
# Snakemake file
#
# 2017.01.13 Added metaspade assembly in addition
#            megahit assembly
# 2017.02.01 Edited to work with Sherlock
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
    qos="bigmem",
    time="16:00:00",
    partition=parameters.ix['subsample_assembly_partition','entry'], 
    mem=parameters.ix['subsample_assembly_memory','entry'],
    kmer=parameters.ix['bulk_megahit_kmerlist','entry']
  threads: int(parameters.ix['subsample_assembly_thread','entry'])
  version: "3.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
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
      shell("""
        echo {scratch}; cd {scratch}; date; pwd; echo
        source activate {python2_env}
        mem=$( echo {params.mem} )000000
        echo 'memory used in Bytes is '$mem 
        megahit_output_dir=megahit_output_{wildcards.subsample}
        source activate {python3_env}
        echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( cat {input_on_scratch[0]| | wc -l ) / 4 ))
        cat {input_on_scratch[2]} {input_on_scratch[3]} > single.fastq
        echo -e 'Starting Assembly: Number of single end reads is: \t'$(( $( wc -l < single.fastq) / 4 ))
        if [ -d $megahit_output_dir ]
        then
        time megahit --k-list {params.kmer} -t {threads} -m $mem --mem-flag 1 --continue -1 {input_on_scratch[0]} {input_on_scratch[1]} -r single.fastq -o $megahit_output_dir --out-prefix megahit
        else
        time megahit --k-list {params.kmer} -t {threads} -m $mem --mem-flag 1 -1 {input_on_scratch[0]} {input_on_scratch[1]} -r single.fastq -o $megahit_output_dir --out-prefix megahit
        fi
        date
        cp $megahit_output_dir/megahit.contigs.fa {output[0]}
        echo 'Assembly Completed'
        source deactivate
        """)
      # Performing Quast
      shell("""
        echo 'Performing metaquast'
        date; pwd; echo
        source activate {python2_env}
        quast_output_dir=metaquast_megahit
        metaquast.py --plots-format svg --gene-finding -m {params.contig_thresh} -t {threads} -o $quast_output_dir {output}
        source deactivate
        """)
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
    qos="normal",
    time="3-0",
    partition=parameters.ix['biosample_assembly_partition','entry'], 
    mem=parameters.ix['biosample_assembly_memory','entry'],
    kmer=parameters.ix['bulk_spades_kmerlist','entry']
  threads: int(parameters.ix['biosample_assembly_thread','entry'])
  version: "3.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
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
      shell("""
        echo {scratch}; cd {scratch}; date; pwd; echo
        source activate {python2_env}
        mem=$( echo {params.mem} | rev | cut -c 4- | rev ) 
        echo 'memory used in Gb is '$mem 
        spades_output_dir=metaSPAdes_output_{wildcards.subsample}
        source activate {python3_env}
        echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( cat {input_on_scratch[0]| | wc -l ) / 4 ))
        if [ -d $spades_output_dir ]
        then
        time metaspades.py -k {params.kmer} -t {threads} --continue -m $mem -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -o $spades_output_dir
        else
        time metaspades.py -k {params.kmer} -t {threads} -m $mem -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -o $spades_output_dir
        fi
        date
        cp $spades_output_dir/contigs.fasta {output}
        echo 'Assembly Completed'
        source deactivate
        """)
      # Performing Quast
      shell("""
        echo 'Performing metaquast'
        date; pwd; echo
        source activate {python2_env}
        quast_output_dir=metaquast_metaSPAdes
        metaquast.py --plots-format svg --gene-finding -m {params.contig_thresh} -t {threads} -o $quast_output_dir {output}
        source deactivate
        """)
      assert(file_empty(output)),"The contig file is empty."
    else:
      print('Input files have size 0')


# This rule has not be updated to work with sherlock yet!!!!!!!!!!
rule align_to_assembly:
  # The output of this rule can be made more informative by including coverage of each contigs.
  input: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq", 
    "{subsample}/contigs.{subsample}.fasta"
  output: "{subsample}/contigCoverage.{subsample}.cnt"
  params: 
    name="bulkBowtie2Align",
    qos="normal",
    time="23:59:59",
    partition=parameters.ix['subsample_bowtie2_partition','entry'], 
    mem=parameters.ix['subsample_bowtie2_memory','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Bowtie2 Alignment of reads back to contigs
    shell("bash {code_dir}/snakehelper_bowtie2align2contig.sh {scratch} {input_on_scratch} \
      {output_on_scratch} {tool_dir} {threads} {wildcards.subsample}")
    cp_from_scratch(output, scratch)


