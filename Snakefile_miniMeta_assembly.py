#################################################################################
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
# 2017.02.21 Updated to be consistent with Sherlock.
#            Removed a lot of helper .sh functions
###################################################################################


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
    qos="normal",
    time="24:00:00",
    partition="normal",
    mem="64000", # don't change this
    contig_thresh=parameters.ix['subsample_contig_thresh','entry']
  threads: 16
  version: "2.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing contig trimming and bowtie2-build
    # The variable work_directory is used in this part
    shell("""
      date; source activate {python2_env}
      cd {scratch};
      cat {input_on_scratch} > combined_contigs.fasta
      python {code_dir}/process_scaffolds.py --lengthThresh {params.contig_thresh} combined_contigs.fasta {output_on_scratch[0]}
      python {code_dir}/process_scaffolds.py --trimHead 150 --trimTail 150 --lengthThresh {params.contig_thresh} {output_on_scratch[0]} trimmed_subsample_contigs.fasta
      source activate {python3_env}
      bowtie2-build --quiet --large-index --threads {threads} -f trimmed_subsample_contigs.fasta subsample_bowtie_{wildcards.id}
      source deactivate; date
      """)
    # Move the rest of the reads back
    cp_from_scratch(output, scratch)

"""
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
    qos="normal",
    time="12:00:00",
    partition=parameters.ix['subsample_bowtie2_partition','entry'], 
    mem=parameters.ix['subsample_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['subsample_contig_thresh','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "2.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing bowtie2 alignment --time --phred33 --un <path> --un-conc <path>
    # bowtie2 can handle fastq.gz files
    shell(""
      date; source activate {python3_env}
      basename=$(echo {input[3]} | cut -d/ -f2 | cut -d. -f1)
      echo $basename
      cd {scratch}
      echo 'Using end-to-end alignment to remove reads'
      bowtie2 --time --phred33 --fast -I 100 -X 2000 -p {threads} --un leftover --un-conc leftover -x $basename -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -U {input_on_scratch[2]} -S alignResults.sam
      mv leftover/un-conc-mate.1 {output_on_scratch[0]}
      mv leftover/un-conc-mate.2 {output_on_scratch[1]}
      mv leftover/un-seqs {output_on_scratch[2]}
      source deactivate
      "")
    # Move the rest of the reads back
    cp_from_scratch(output, scratch)
"""


rule minimeta_spade_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2017.01.12 input is edited to use YAML file and combined subsample contigs
  input: 
    expand("{subsample}/leftover_reads_P1.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/leftover_reads_P2.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/leftover_reads_S.{subsample}.fastq", subsample=subsampleIDs),
    expand("Combined_Analysis/subsample_contigs.{id}.fasta", id=biosample)
  output: # metaquast does not require output file, just output dir
    "miniMetaAssembly_reads.yaml", # in the root work_directory because of spades paths
    "Combined_Analysis/minimeta_contigs.{id}.fasta", 
    "Combined_Analysis/minimeta_scaffolds.{id}.fasta" 
  params:    
    name="minimeta_spade_assembly",
    qos="normal",
    time="7-0",  # means 7 days
    partition=parameters.ix['biosample_assembly_partition','entry'], 
    mem=parameters.ix['biosample_assembly_memory','entry'],
    contig_thresh=parameters.ix['biosample_contig_thresh','entry'],
    kmer=parameters.ix['biosample_spades_kmerlist','entry'] # "21,33,44,77,99"
  threads: int(parameters.ix['biosample_assembly_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    print(work_directory)
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
    print('Completed creating YAML file')
    # Performing spade assembly, not using --meta
    # quast.py --plots-format svg --gene-finding -meta -m {params.contig_thresh} -t {threads} -o {scratch}/quast_output {output[1]}
    shell("""
      echo {scratch}; date; pwd; echo
      source activate {python2_env}
      mem=$( echo {params.mem} | rev | cut -c 4- | rev )
      echo 'memory used in Gb is '$mem
      spades_output_dir=Combined_Analysis/spades_minimetaAssembly_{wildcards.id}
      source activate {python3_env}
      if [ -d $spades_output_dir ]
      then
      time spades.py --only-assembler -k {params.kmer} -t {threads} --continue --sc -m $mem --dataset {output[0]} -o $spades_output_dir
      else
      time spades.py --only-assembler -k {params.kmer} -t {threads} --sc -m $mem --dataset {output[0]} -o $spades_output_dir
      fi
      cp $spades_output_dir/contigs.fasta {output[1]}
      cp $spades_output_dir/scaffolds.fasta {output[2]}
      source activate {python2_env}
      metaquast.py --plots-format svg --gene-finding -m {params.contig_thresh} -t {threads} -o $spades_output_dir/metaquast_output {output[1]}
      date; source deactivate
      """)
    assert(file_empty(output)),"Either the contig or scaffold file is empty."


