##############################################################
# Snakemake rules associated with assembly
# blast, quast of a single subsample.
# this file must be included into another 
# Snakemake file
#
# 2017.01.31 Adapted this file for sherlock.stanford.edu
#            scratch = os.environ["LOCAL_SCRATCH"]
#            each "normal" partition core contains 4GB memory
##############################################################

# rules

rule subsample_spade_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2015.11.30 This rule is still using splitNcluster outputs.
  #            These need to be changed to the P1 P2 S1 S2 from trim
  # 2015.12.01 Changes were made to use trimmed reads directly
  input: 
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
    qos="normal",
    time="12:00:00",
    partition=parameters.ix['subsample_assembly_partition','entry'], 
    mem=parameters.ix['subsample_assembly_memory','entry'],
    kmer=parameters.ix['subsample_spades_kmerlist','entry'] # "21,33,55,77,99"
  threads: int(parameters.ix['subsample_assembly_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
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
      # Assembly using spades V3.9. Only include output_on_scratch[0:2] here 
      # because the snakemake_spadeassembly file is not written to handle the 
      # other two outputs. However you can copy those back later from SCRATCH
      # by the way, [0:2] actually only contains 0th and 1st entry
      # 2015.12.01 Edited to use trimmed reads directly
      # 2017.01.31 Edited so that kmer is now a variable string
      shell("""
        echo {scratch}; cd {scratch}
        source activate {python3_env}
        spadesMemory=$( echo {params.mem} | rev | cut -c 4- | rev )
        echo 'memory used in GB is: '$spadesMemory
        echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( wc -l < {input_on_scratch[0]}) / 4 ))
        cat {input_on_scratch[2]} {input_on_scratch[3]} > single.fastq
        echo -e 'Starting Assembly: Number of single end reads is: \t'$(( $( wc -l < single.fastq ) / 4 ))
        subsample_spades_dir=spades_output_{wildcards.subsample}
        time spades.py -k {params.kmer} -t {threads} --sc --careful -m $spadesMemory -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -s single.fastq -o $subsample_spades_dir
        cp $subsample_spades_dir/contigs.fasta {output_on_scratch[0]}
        cp $subsample_spades_dir/scaffolds.fasta {output_on_scratch[1]}
        mv $subsample_spades_dir/corrected/*1.*.00.0_0.cor.fastq.gz {output_on_scratch[2]}
        mv $subsample_spades_dir/corrected/*2.*.00.0_0.cor.fastq.gz {output_on_scratch[3]}
        mv $subsample_spades_dir/corrected/*unpaired.00.0_0.cor.fastq.gz {output_on_scratch[4]}
        mv $subsample_spades_dir/spades.log {output_on_scratch[5]}
        mv $subsample_spades_dir/params.txt {output_on_scratch[6]}
        echo; ls; date
        # source deactivate
        """)
      cp_from_scratch(output[2:], scratch)
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
      print('Assembly Completed')
    else:
      print('Input files have size 0')


rule subsample_remove_reads:
  # The output of this rule can be made more informative by including contig cov.
  # 2017.02.01 Changed to use corrected reads for alignment, bowtie2 handles .gz
  input: 
    # "{subsample}/P1.{subsample}.fastq", 
    # "{subsample}/P2.{subsample}.fastq",
    "{subsample}/P1_corrected.{subsample}.fastq.gz", 
    "{subsample}/P2_corrected.{subsample}.fastq.gz",
    "{subsample}/S_corrected.{subsample}.fastq.gz",
    "{subsample}/contigs.{subsample}.fasta"
  output:
    temp("{subsample}/leftover_reads_P1.{subsample}.fastq"),
    temp("{subsample}/leftover_reads_P2.{subsample}.fastq"),
    temp("{subsample}/leftover_reads_S.{subsample}.fastq")
  params: 
    name="subsample_remove_reads",
    qos="normal",
    time="2:00:00",
    partition=parameters.ix['subsample_bowtie2_partition','entry'], 
    mem=parameters.ix['subsample_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['subsample_contig_thresh','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Bowtie2 Alignment of reads back to contigs
    shell("""
      echo {scratch}; cd {scratch}; date
      source activate {python2_env}
      python {code_dir}/process_scaffolds.py --lengthThresh {params.contig_thresh} {input_on_scratch[3]} trimmed_subsample_contigs.fasta
      source activate {python3_env}
      echo 'Building bowtie2 indices'
      bowtie2-build --quiet --threads {threads} -f trimmed_subsample_contigs.fasta subsampleContigs
      date; echo; du -sh; echo; echo 'Aligning reads using bowtie2'
      bowtie2 --time --phred33 --very-sensitive -I 100 -X 2000 -p {threads} --un leftoverReadsSingle --un-conc leftoverReads -x subsampleContigs -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -U {input_on_scratch[2]} -S alignedResults.sam
      mv leftoverReads.1 {output_on_scratch[0]}
      mv leftoverReads.2 {output_on_scratch[1]}
      mv leftoverReadsSingle {output_on_scratch[2]}
      echo; date; ls; echo; du -sh; echo; source deactivate
      """)
    cp_from_scratch(output, scratch)


rule subsample_quast:
  input: "{folder}/contigs.{k}.fasta"
  output: "{folder}/quast_report.{k}.txt"
  params: 
    name="subsample_quast",
    qos="normal",
    time="30:00",
    partition="normal,quake,owners", 
    mem="8000",
    contig_thresh=int(parameters.ix['biosample_contig_thresh','entry'])
  threads: 2
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing Quast v4.4 with biosample_contig_thresh value
    # quast is in the anaconda2 env, not quast.py just quast
    # could use options --plots-format svg --gene-finding --meta
    shell("""
      echo {scratch}; cd {scratch}
      source activate {python2_env}
      quast -m {params.contig_thresh} -t {threads} -o subsample_quast_output {input_on_scratch}
      mv subsample_quast_output/report.txt {output_on_scratch}
      source deactivate
      """)
    cp_from_scratch(output, scratch)

    
