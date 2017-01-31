###############################################
# Snakemake rules associated with analyzing 
# super contigs in conjunction with subsamples
# this file must be included into another 
# Snakemake file 
###############################################   

def merge_superContig_alignment_results(input, output):
  """
  input: a list of files in the form of subsample_superContig_alignmentReport.ILxxxxxx.txt
  output: table containing total coverage average coverage of each contig in each sample
          normalized by contig length
  """
  num_input = len(input) - 1
  assert(num_input > 0),"Input has no valid files."
  # The first input is the super contig file. 
  # Each super_contig must contain exactly 2 lines, The first line must begin with '>' 
  # and the second line must contain only the sequence and the entire sequence.
  # Also all contig names in the input file must be unique
  superContig_file = input[0]
  with open(superContig_file,'r') as f:
    contig_name = []
    contig_length = []
    for l in f:
      if l[0] == '>':
        assert(len(contig_name) == len(contig_length)),"Missed a contig length."
        contig_name.append(l.split()[0][1:]) # ignores the first '>'
      else:
        assert(len(contig_name) == len(contig_length)+1),"Missed a contig name."
        contig_length.append(len(l.split()[0]))
  assert(len(contig_name)==len(contig_length)),"Contig name array and contig length array have different number of elements."
  # Turn contig length into a dictionary and coverage into a different dictionary
  superContig_length = {contig_name[i] : contig_length[i] for i in range(len(contig_length))}
  superContig_coverage = {contig_name[i] : [0 for j in range(num_input)] for i in range(len(contig_length))}
  # Each input file contains 2 columns. The first column is the supercontig name
  # the second column is the number of basepairs covered in total
  filenames = input[1:]
  subsamples = []
  for ind in range(len(filenames)):
    filename = filenames[ind]
    print("Currently processing file: "+filename)
    # file name should be in the form of Combined_Analysis/subsample_superContig_alignmentReport.ILxxxxx.txt
    print("The corresponding sample name is: "+filename.split('.')[1])
    subsamples.append(filename.split('.')[1])
    with open(filename,'r') as f:
      for l in f:
        # this is a hack because all the subsample alignment reports have an empty line at the beginning.
        if l.split(): # if the line is not an empty line
          l = l.split()
          superContig_coverage[l[0]][ind] = float(l[1])/float(superContig_length[l[0]])
  # Processing output file. output should be a list of only one string
  print(subsamples)
  with open(output[0],'w') as f:
    subsamples.insert(0, '')
    t = f.write('\t'.join(subsamples) + '\n')
    for key in superContig_coverage:
      # add the supercontig name in front of all the coverages
      # print(key + '\t' + '\t'.join([str(x) for x in superContig_coverage[key]]))
      # separate by tab and add new line at the end (need to turn int list to str list first)
      t = f.write(key + '\t' + '\t'.join([str(x) for x in superContig_coverage[key]]) + '\n')

"""
rule subsampleReads_Align2SuperContigs:
  input: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq", 
    expand("Combined_Analysis/super_contigs.{id}.fasta", id=biosample)
  output: 
    temp("Combined_Analysis/subsample_superContig_alignmentReport.{subsample}.txt"),
    "Combined_Analysis/subsample_superContig_alignment.{subsample}.bam",
    "Combined_Analysis/subsample_superContig_alignment.{subsample}.bcf",
    "Combined_Analysis/subsample_superContig_alignment.{subsample}.vcf"
  params:
    name="subsampleReads_Align2SuperContigs",
    partition=parameters.ix['subsample_bowtie2_partition','entry'],
    mem=parameters.ix['subsample_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location                                                        
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Bowtie2 Alignment of reads back to super_contigs (contigs with new names)
    shell("bash {code_dir}/snakehelper_bowtie2SampleSupercontig.sh {scratch} {input_on_scratch} \
      {output_on_scratch} {tool_dir} {code_dir} {params.contig_thresh} {threads}")
    cp_from_scratch(output, scratch)
"""

rule merge_superContig_Alignment:
  input: 
    "Combined_Analysis/super_contigs.{id}.fasta",
    expand("Combined_Analysis/subsample_superContig_alignmentReport.{subsample}.txt", subsample=subsampleIDs)
  output: "Combined_Analysis/super_contigs.{id}.alignment_report.txt"
  params:
    name="merge_superContig_Alignment",
    partition="general",
    mem="5000", # don't change
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: 1
  version: "1.0"
  run: 
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    contig_on_scratch = names_on_scratch(["super_contig_subset.fasta"], scratch)
    cp_to_scratch(input, scratch)
    # Perform organization of contigs    
    shell("set +u; source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/ &&\
      python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[0]} {contig_on_scratch[0]}")
    input_on_scratch[0] = contig_on_scratch[0]
    merge_superContig_alignment_results(input_on_scratch, output_on_scratch)
    cp_from_scratch(output, scratch)


