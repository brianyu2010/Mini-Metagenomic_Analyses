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


def addToContigName(inputFile, front_flag, string2add, outputFile):
  """
  inputFile: the file name containing original contigs
  front_flag: if true, add string to front ie. >string_xxx else
              add string to the end
  string2add: self explanatory, does not included spacer '_'
  outputFile: name of the output contig file
  description:  
  This function also turns all the spaces in megahit contig names to '_'
  """
  with open(inputFile,'r') as f1, open(outputFile,'w') as f2:
    fasta_label = []
    fasta_seq = []
    counter = 0
    for l in f1:
      # turn all spaces in the line into '_' (for megahit contigs)
      if ' ' in l:
        l = l.replace(' ','_')
      # check if the line begins with '>', meaning it's a header
      if '>' in l:
        counter = counter + 1
        # check if there is no '\n'
        if '\n' not in l:
          if front_flag:
            l = '>' + string2add + '_' + l[1:] + '\n'
          else:
            l = l + '_' + string2add + '\n'
        # if there is a '\n' character
        else:
          if front_flag:
            l = '>' + string2add + '_' + l[1:]
          else:
            l = l.split()[0] + '_' + string2add + '\n'
        # Append to label
        fasta_label.append(l)
      else: # not a label then it's a sequence
        if '\n' not in l:
          l.append('\n')
        # if length of names > length of seq then add l to seq
        # because this is the first line of sequence
        if len(fasta_label) > len(fasta_seq):
          fasta_seq.append(l)
        else:
          fasta_seq[-1] = fasta_seq[-1][0:-1] + l
    assert(len(fasta_seq) == len(fasta_label))
    print('Number of lines is %d' % (counter))
    # Write to file
    for i in range(len(fasta_label)):
      t = f2.write('%s%s' %(fasta_label[i], fasta_seq[i]))
  print('Updated contig names in file ' + inputFile + ' and saved in ' + outputFile)


rule readsAndContigs_spades_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2017.01.12 input is edited to use YAML file and combined subsample contigs
  # 2017.02.01 combined YAML file creation, assembly, quast into one rule
  input: 
    expand("{subsample}/P1_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/P2_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/S_corrected.{subsample}.fastq.gz", subsample=subsampleIDs),
    expand("{subsample}/contigs.{subsample}.fasta", subsample=subsampleIDs)
  output: # for metaquast, no output is required; all output need the same wildcards
    "corrected_readsAndContigs.{id}.yaml",
    "Combined_Analysis/readsAndContigsAssembly_contigs.{id}.fasta", 
    "Combined_Analysis/readsAndContigsAssembly_scaffolds.{id}.fasta",
    "Combined_Analysis/quast_report_readsAndContigs.{id}.txt"
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
      spades_output_dir=Combined_Analysis/spades_readsAndContigsAssembly_{wildcards.id}
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
      # metaquast.py --plots-format svg --gene-finding -m {params.contig_thresh} -t {threads} -o $spades_output_dir/metaquast_output {output[1]}
      quast.py -m {params.contig_thresh} -t {threads} -o $spades_output_dir/quast_output {output[1]}
      cp $spades_output_dir/quast_output/report.txt {output[3]}
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
    "corrected_readsOnly.{id}.yaml",
    "Combined_Analysis/readsOnlyAssembly_contigs.{id}.fasta", 
    "Combined_Analysis/readsOnlyAssembly_scaffolds.{id}.fasta",
    "Combined_Analysis/quast_report_readsOnly.{id}.txt"
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
      spades_output_dir=Combined_Analysis/spades_readsOnlyAssembly_{wildcards.id}
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
      # metaquast.py --plots-format svg --gene-finding -m {params.contig_thresh} -t {threads} -o $spades_output_dir/metaquast_output {output[1]}
      quast.py -m {params.contig_thresh} -t {threads} -o $spades_output_dir/quast_output {output[1]}
      cp $spades_output_dir/quast_output/report.txt {output[3]}
      date; source deactivate
      """)
    assert(file_empty(output)),"Either the contig or scaffold file is empty."


if bulk_flag =='Yes' or bulk_flag == 'yes' or bulk_flag == 'Y' or bulk_flag == 'y':
  # rename each super contig with super contig code
  rule make_superContigs_bulk:
    input:
      expand("{bulksample}/megahit_contigs.{bulksample}.fasta", bulksample=bulksampleIDs),
      expand("{bulksample}/metaSPAdes_contigs.{bulksample}.fasta", bulksample=bulksampleIDs),
      "{folder}/minimeta_contigs.{id}.fasta"
    output:
      "{folder}/super_contigs.{id}.fasta"
    params:
      name="make_super_contigs_bulk",
      qos="normal",
      time="16:00:00",
      partition="normal,quake",
      mem="126000",
      contig_thresh=parameters.ix['biosample_contig_thresh','entry']
    threads: 20
    version: "2.0"
    run:
      # Manage files and obtain scratch location
      scratch = os.environ["LOCAL_SCRATCH"]
      input_on_scratch = names_on_scratch(input, scratch)
      # The set +u; prevents unbound variable errors 2016.08.18
      # down samples all contigs
      for i in range(len(input)):
        inputFileName = input[i]
        scratchFileName = input_on_scratch[i]
        shell("""
          source activate {python2_env}
          python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {inputFileName} {scratch}/temp.fasta
          source deactivate
          echo {inputFileName}
          echo {scratchFileName}
          """)
        # Add contig names to front
        addToContigName(scratch+'/temp.fasta', 1, 'SuperContig', scratchFileName)
      # Combine contigs
      shell("""
        cat {input_on_scratch} > {output}
        echo 'SuperContig construction completed'
        """)
      # Making metaquast and quast for combined supercontigs ONLY when bulk is included
      shell("""
        echo; date; pwd 
        source activate {python2_env}
        metaquast_output_dir={wildcards.folder}/metaquast_superContig
        quast_output_dir={wildcards.folder}/quast_superContig
        # --gene-finding --meta has expired
        metaquast.py --plots-format svg --max-ref-number 200 -m {params.contig_thresh} -t {threads} -o $metaquast_output_dir {output}
        # --gene-finding --meta has expired. No longer do gene finding
        quast.py -m {params.contig_thresh} -t {threads} -o $quast_output_dir {output}
        date; echo; source deactivate
        """)
      assert(file_empty(output)),"Supercontig file is empty."
else:
  # rename each super contig with super contig code
  rule make_superContigs_miniMetaOnly:
    input:
      "{folder}/minimeta_contigs.{id}.fasta"
    output:
      "{folder}/super_contigs.{id}.fasta"
      # "{folder}/super_contigs_name.{id}.txt" # no longer needed
    params:
      name="make_super_contigs_miniMetaOnly",
      qos="normal",
      time="2:00:00",
      partition="normal,quake",
      mem="50000",
      contig_thresh=parameters.ix['biosample_contig_thresh','entry']
    threads: 4
    version: "2.0"
    run:
      # Manage files and obtain scratch locations
      scratch = os.environ["LOCAL_SCRATCH"]
      input_on_scratch = names_on_scratch(input, scratch)
      # down samples all contigs, and move to scratch
      for i in range(len(input)):
        inputFileName = input[i]
        scratchFileName = input_on_scratch[i]
        shell("""
          source activate {python2_env}
          python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {inputFileName} {scratch}/temp.fasta
          source deactivate
          echo {inputFileName}
          echo {scratchFileName} 
          """)
        # Add contig names to front
        addToContigName(scratch+'/temp.fasta', 1, 'SuperContig', scratchFileName)
      # Combine contigs
      shell("""
        cat {input_on_scratch} > {output}
        echo 'SuperContig construction completed'
        """)
      assert(file_empty(output)),"Supercontig file is empty."

