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


###############################################
# Snakemake rules associated with analyzing 
# super contigs in conjunction with subsamples
# to compute size of genome from each subsample
# this file must be included into another 
# Snakemake file 
#
# 2016.08.30 Created Brian Yu
###############################################   

rule subsampleSuperContigPileup:
  input: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq", 
    expand("Combined_Analysis/super_contigs.{id}.fasta", id=biosample)
  output: 
    temp("Combined_Analysis/subsample_superContig_mpileupReport.{subsample}.txt"),
    temp("Combined_Analysis/subsampleAlign2Supercontig.{subsample}.bam"),
    temp("Combined_Analysis/subsample_superContig_alignmentReport.{subsample}.txt")
  params:
    name="subsampleSuperContigPileup",
    partition="general", #parameters.ix['subsample_bowtie2_partition','entry'],
    mem="21000", #parameters.ix['subsample_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: 4 #int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Bowtie2 Alignment of reads back to super_contigs (contigs with new names)
    # This is perhaps not the best way, use the line above instead
    # awk 'a1==$1 {a2+=$4; next} {print a1,"\t", a2; a1=$1; a2=$4} END {print a1,"\t",a2}' alignResults.pile > {output_on_scratch[2]}  
    shell("""cp {code_dir}/threshold_scaffolds.py {scratch}
      cd {scratch}; set +u
      source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
      python threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[2]} temp.fasta
      {tool_dir}/bowtie2-2.2.6/bowtie2-build -f temp.fasta spadeContigs
      {tool_dir}/bowtie2-2.2.6/bowtie2 --very-sensitive-local -I 100 -X 2000 -p {threads} -t -x spadeContigs -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -S alignResults.sam
      {tool_dir}/samtools-1.3/samtools view -b -o alignResults.bam alignResults.sam
      {tool_dir}/samtools-1.3/samtools sort -o alignResults_sorted.bam alignResults.bam
      cp alignResults_sorted.bam {output_on_scratch[1]}
      {tool_dir}/samtools-1.3/samtools index alignResults_sorted.bam
      {tool_dir}/samtools-1.3/samtools mpileup -f temp.fasta -o alignResults.pile alignResults_sorted.bam 
      awk '$4>=5 {{print}}' alignResults.pile | cut -f 1 | sort | uniq -c | awk '{{print $2,"\t",$1}}' > {output_on_scratch[0]}
      echo nextLineIsTheProblem
      awk '$4>=1 {{print}}' alignResults.pile | cut -f 1 | sort | uniq -c | awk '{{print $2,"\t",$1}}' > {output_on_scratch[2]}
      echo Process_Completed
      """)
    cp_from_scratch(output, scratch)



rule VCFgenerate:
  input:
    "Combined_Analysis/super_contigs.{id}.fasta",
    expand("Combined_Analysis/subsampleAlign2Supercontig.{subsample}.bam", subsample=subsampleIDs)
  output:
    "Combined_Analysis/subsample_variants.{id}.vcf",
    "Combined_Analysis/subsample_supercontigCoverage.{id}.pile"
  params:
    name="VCFgenerate",
    partition="long",
    mem="60000",
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: 2
  version: "1.0"
  run: 
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    contig_on_scratch = names_on_scratch(["super_contig_subset.fasta"], scratch)
    bamfiles_on_scratch = input_on_scratch[1:]
    cp_to_scratch(input, scratch)
    # Perform organization of contigs    
    shell("""cd {scratch}; set +u; source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
      python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[0]} {contig_on_scratch[0]}
      {tool_dir}/samtools-1.3/samtools mpileup -g -t DP -f {contig_on_scratch[0]} -o subsample_alignment.bcf {bamfiles_on_scratch}
      echo BCF_file_generation_completed
      {tool_dir}/samtools-1.3/samtools mpileup -f {contig_on_scratch[0]} -o {output_on_scratch[1]} {bamfiles_on_scratch}
      {tool_dir}/bcftools-1.3/bcftools call -c -v -o {output_on_scratch[0]} subsample_alignment.bcf 
      echo VCF_file_generation_completed
      """)
    cp_from_scratch(output, scratch)
  


rule subsampleGenomeSize:
  input: 
    "Combined_Analysis/super_contigs.{id}.fasta",
    expand("Combined_Analysis/subsample_superContig_mpileupReport.{subsample}.txt", subsample=subsampleIDs)
  output: 
    "Combined_Analysis/super_contigs.{id}.subsampleGenomeSize.txt"
  params:
    name="subsampleGenomeSize",
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
    shell("""set +u; source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
      python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[0]} {contig_on_scratch[0]}""")
    input_on_scratch[0] = contig_on_scratch[0]
    merge_superContig_alignment_results(input_on_scratch, output_on_scratch)
    cp_from_scratch(output, scratch)

