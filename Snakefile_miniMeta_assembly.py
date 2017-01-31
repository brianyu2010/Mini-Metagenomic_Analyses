###############################################
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
# 2017.01.21 Modified biosample assembly and created this
#            miniMeta_assembly.py. It does these steps:
#            1. Combine all subsample contigs, threshold
#            2. Remove 200 bp from each end and align reads from each sub-sample
#            3. remove reads that align and assemble with --trusted-contigs
###############################################


rule combine_threshold_subsample_contigs:
  input:
    expand("{subsample}/contigs.{subsample}.fasta", subsample=subsampleIDs)
  output:
    "Combined_Analysis/subsample_contigs.{id}.fasta",
    "Combined_Analysis/subsample_bowtie_{id}.1.bt2l",
    "Combined_Analysis/subsample_bowtie_{id}.2.bt2l",
    "Combined_Analysis/subsample_bowtie_{id}.3.bt2l",
    "Combined_Analysis/subsample_bowtie_{id}.4.bt2l",
    "Combined_Analysis/subsample_bowtie_{id}.rev.1.bt2l",
    "Combined_Analysis/subsample_bowtie_{id}.rev.2.bt2l"
  params:
    name="combine_threshold_subsample_contigs",
    partition="long",
    mem="63000", # don't change this
    contig_thresh=parameters.ix['subsample_contig_thresh','entry']
  threads: 12
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing contig trimming and bowtie2-build
    # The variable work_directory is used in this part
    shell("""
      date
      cp {code_dir}/threshold_scaffolds.py {scratch}
      cp {code_dir}/process_scaffolds.py {scratch}
      cd {scratch}; set +u
      cat {input_on_scratch} > combined_contigs.fasta
      source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
      python process_scaffolds.py --lengthThresh {params.contig_thresh} combined_contigs.fasta {output_on_scratch[0]}
      python process_scaffolds.py --trimHead 150 --trimTail 150 --lengthThresh {params.contig_thresh} {output_on_scratch[0]} trimmed_subsample_contigs.fasta
      {tool_dir}/bowtie2-2.2.6/bowtie2-build --large-index -f trimmed_subsample_contigs.fasta subsample_bowtie_{wildcards.id}
      date
      """)
    # Move the rest of the reads back
    cp_from_scratch(output, scratch)


# This rule also needs pre-built bowtie index, not part of input
rule remove_aligned_reads:
  input:
    "{subsample}/P1_corrected.{subsample}.fastq.gz",
    "{subsample}/P2_corrected.{subsample}.fastq.gz",
    "{subsample}/S_corrected.{subsample}.fastq.gz",
    expand("Combined_Analysis/subsample_bowtie_{id}.1.bt2l", id=biosample),
    expand("Combined_Analysis/subsample_bowtie_{id}.2.bt2l", id=biosample),
    expand("Combined_Analysis/subsample_bowtie_{id}.3.bt2l", id=biosample),
    expand("Combined_Analysis/subsample_bowtie_{id}.4.bt2l", id=biosample),
    expand("Combined_Analysis/subsample_bowtie_{id}.rev.1.bt2l", id=biosample),
    expand("Combined_Analysis/subsample_bowtie_{id}.rev.2.bt2l", id=biosample)
  output:
    "{subsample}/leftover_reads_P1.{subsample}.fastq",
    "{subsample}/leftover_reads_P2.{subsample}.fastq",
    "{subsample}/leftover_reads_S.{subsample}.fastq"
  params:
    name="remove_aligned_reads",
    partition=parameters.ix['subsample_bowtie2_partition','entry'], 
    #partition="general",
    mem=parameters.ix['subsample_bowtie2_memory','entry'],
    #mem="10500"
    contig_thresh=parameters.ix['subsample_contig_thresh','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  #threads: 2
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing bowtie2 alignment
    # Trim contigs by 150 bp at both ends first
    # align, extract unaligned reads, turn bam to fastq
    # bowtie2 can handle fastq.gz files
    shell("""
      date
      basename=$(echo {input[3]} | cut -d/ -f2 | cut -d. -f1)
      echo $basename
      cd {scratch}
      {tool_dir}/bowtie2-2.2.6/bowtie2 --fast-local -I 100 -X 2000 -p {threads} -t -x $basename -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -U {input_on_scratch[2]} -S alignResults.sam
      {tool_dir}/samtools-1.3/samtools view -b -F 0x002 -o unmapped.bam alignResults.sam
      {tool_dir}/samtools-1.3/samtools fastq -1 {output_on_scratch[0]} -2 {output_on_scratch[1]} unmapped.bam > {output_on_scratch[2]}
      """)
    # Move the rest of the reads back
    cp_from_scratch(output, scratch)



rule create_miniMeta_assembly_YAML:
  input:
    expand("{subsample}/leftover_reads_P1.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/leftover_reads_P2.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/leftover_reads_S.{subsample}.fastq", subsample=subsampleIDs),
    expand("Combined_Analysis/subsample_contigs.{id}.fasta", id=biosample)
  output:
    # This is in the root work_directory
    "miniMetaAssembly_reads.yaml"
  params:
    name="create_miniMeta_assembly_YAML",
    partition="general",
    mem="1000" # don't change this
  threads: 1
  version: "1.0"
  run:
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
      # Write trusted contigs for this assembly, there should only be one file
      d = f.write('  {\n')
      d = f.write('    type: trusted-contigs,\n')
      d = f.write('    single reads: [\n')
      d = f.write('      ' + input[-1] + '\n')
      d = f.write('    ]\n')
      d = f.write('  }\n')
      d = f.write(']')



rule minimeta_spade_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2017.01.12 input is edited to use YAML file and combined subsample contigs
  input: 
    "miniMetaAssembly_reads.yaml", # in the root work_directory because of spades paths
    expand("{subsample}/leftover_reads_P1.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/leftover_reads_P2.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/leftover_reads_S.{subsample}.fastq", subsample=subsampleIDs),
    expand("Combined_Analysis/subsample_contigs.{id}.fasta", id=biosample)
  output: 
    "Combined_Analysis/minimeta_contigs.{id}.fasta", 
    "Combined_Analysis/minimeta_scaffolds.{id}.fasta",
    "Combined_Analysis/quast_report_miniMeta.{id}.txt"
  params:    
    name="minimeta_spade_assembly", 
    partition=parameters.ix['biosample_assembly_partition','entry'], 
    mem=parameters.ix['biosample_assembly_memory','entry']
  threads: int(parameters.ix['biosample_assembly_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    print(work_directory)
    # Assembly using spades V3.5.0 much better kmer normalization
    # shell("bash {code_dir}/snakehelper_bigmemAssembly.sh {scratch} {input} {tool_dir} {params.mem} {threads} {work_directory}/Combined_Analysis/spade_output_{wildcards.id}")
    shell("bash {code_dir}/snakehelper_bigmemYAMLAssembly.sh {scratch} {input[0]} {tool_dir} {params.mem} {threads} {work_directory}/Combined_Analysis/spades_minimetaAssembly_{wildcards.id}")
    shell("cp Combined_Analysis/spades_minimetaAssembly_{wildcards.id}/contigs.fasta {output[0]}")
    shell("cp Combined_Analysis/spades_minimetaAssembly_{wildcards.id}/scaffolds.fasta {output[1]}")
    shell("""
      set +u
      source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
      python {tool_dir}/quast-3.2/quast.py -o {scratch}/quast_output {output[0]}
      cp {scratch}/quast_output/report.txt {output[2]}
      """)
    assert(file_empty(output)),"Either the contig or scaffold file is empty."


