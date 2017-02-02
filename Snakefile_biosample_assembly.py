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
#            minimeta contigs.
#            Made it consistent with sherlock
##################################################################
      
rule readsAndContigs_spades_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2017.01.12 input is edited to use YAML file and combined subsample contigs
  # 2017.02.01 combined YAML file creation, assembly, quast into one rule
  input: 
    expand("{subsample}/P1_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/P2_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/S_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/contigs.{subsample}.fasta", subsample=subsampleIDs)
  output: # for metaquast, no output is required
    "corrected_readsAndContigs.yaml" # in the root work_directory because of spades paths
    "Combined_Analysis/readsAndContigsAssembly_contigs.{id}.fasta", 
    "Combined_Analysis/readsAndContigsAssembly_scaffolds.{id}.fasta"
  params:    
    name="readsAndContigs_spades_assembly",
    qos="normal",
    time="7-0", # means 7 days
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
      # Write untrusted contigs containing all subsample contigs
      d = f.write('  {\n')
      d = f.write('    type: untrusted-contigs,\n')
      d = f.write('    single reads: [\n')
      d = f.write('      ' + ',\n      '.join(input[(3*numfiles):(4*numfiles)]) + '\n')
      d = f.write('    ]\n')
      d = f.write('  }\n')
      d = f.write(']')
    print('Completed creating YAML file')
    # Performing Assembly
    shell("""
      echo {scratch}; date; pwd; echo
      source activate {python2_env}
      mem=$( echo {params.mem} | rev | cut -c 4- | rev )
      echo 'memory used in Gb is '$mem 
      spades_output_dir=Combined_Analysis/spade_readsAndContigsAssembly_{wildcards.id}
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


rule readsOnly_spades_assembly:
  input: 
    expand("{subsample}/P1_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/P2_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/S_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/contigs.{subsample}.fasta", subsample=subsampleIDs) # not used
  output:
    "corrected_readsOnly.yaml" # in the root work_directory because of spades paths
    "Combined_Analysis/readsOnlyAssembly_contigs.{id}.fasta", 
    "Combined_Analysis/readsOnlyAssembly_scaffolds.{id}.fasta"
  params:    
    name="readsOnly_spades_assembly",
    qos="normal",
    time="7-0", # means 7 days
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
      d = f.write('  }\n')
      d = f.write(']')
    print('Completed creating YAML file')
    # Performing Assembly
    shell("""
      echo {scratch}; date; pwd; echo
      source activate {python2_env}
      mem=$( echo {params.mem} | rev | cut -c 4- | rev )
      echo 'memory used in Gb is '$mem 
      spades_output_dir=Combined_Analysis/spade_readsOnlyAssembly_{wildcards.id}
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
    

# rename each super contig with super contig code
rule make_superContigs:
  input:
    # "{folder}/combined_contigs.{id}.fasta"
    "{folder}/minimeta_contigs.{id}.fasta"
  output:
    "{folder}/super_contigs.{id}.fasta",
    "{folder}/super_contigs_name.{id}.txt"
  params:
    name="make_super_contigs",
    qos="normal",
    time="30:00",
    partition="normal",
    mem="4000"
  threads: 1
  version: "2.0"
  run: 
    # The set +u; prevents unbound variable errors 2016.08.18
    shell(""" 
      python --version; 
      source activate {python2_env}; 
      python --version
      python {code_dir}/snakehelper_construct_superContigs.py {input} -o {output[0]} -l {output[1]}
      source deactivate
      echo 'SuperContig construction completed'
      """)
    assert(file_empty(output)),"Either the supercontig or supercontig names are empty."


# The blast database has not been setup yet so don't use this rule yet
"""
rule biosample_BLAST:
  input: "Combined_Analysis/super_contigs.{id}.fasta"
  output: "Combined_Analysis/BlastResults.{id}.txt"
  params:
    name="biosample_blast_results",
    partition="normal",
    mem="64000",
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: 16
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Blast contigs
    assert(file_empty([input[0]])),"Input Contig file is empty."
    shell(""
      echo {scratch}
      source activate {python3_env}


      bash {code_dir}/snakehelper_blast.sh {scratch} {input_on_scratch} {output_on_scratch} {tool_dir} {code_dir} {threads} {params.contig_thresh}
      source deactivate
      "")
    assert(file_empty([output_on_scratch[0]])),"Blast output file is empty."
    cp_from_scratch(output, scratch)

"""
