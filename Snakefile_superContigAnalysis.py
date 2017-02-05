#######################################################
# Snakemake rules associated with analyzing 
# super contigs in conjunction with subsamples
# this file must be included into another 
# Snakemake file 
#
# 2017.02.01 Updated to work on the sherlock cluster
#######################################################

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


if bulk_flag =='Yes' or bulk_flag == 'yes' or bulk_flag == 'Y' or bulk_flag == 'y':
  rule merge_superContig_alignment_withBulk:
    input: 
      "Combined_Analysis/super_contigs.{id}.fasta",
      expand("Combined_Analysis/subsample_superContig_alignmentReport.{subsample}.txt", subsample=subsampleIDs),
      expand("Combined_Analysis/subsample_superContig_alignmentReport.{bulksample}.txt", bulksample=bulksampleIDs)
    output: "Combined_Analysis/super_contigs.{id}.alignment_report.txt"
    params:
      name="merge_superContig_alignment_withBulk",
      qos="normal",
      time="2:00:00",
      partition="normal",
      mem="4000", # don't change
      contig_thresh=parameters.ix['biosample_contig_thresh','entry']
    threads: 1
    version: "2.0"
    run: 
      # Managing files and obtain scratch location
      scratch = os.environ["LOCAL_SCRATCH"]
      input_on_scratch = names_on_scratch(input, scratch)
      output_on_scratch = names_on_scratch(output, scratch)
      contig_on_scratch = names_on_scratch(["super_contig_subset.fasta"], scratch)
      cp_to_scratch(input, scratch)
      # Perform organization of contigs    
      shell("""
        source activate {python2_env}
        python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[0]} {contig_on_scratch[0]}
        """)
      input_on_scratch[0] = contig_on_scratch[0]
      merge_superContig_alignment_results(input_on_scratch, output_on_scratch)
      cp_from_scratch(output, scratch)
else:
  rule merge_superContig_alignment_miniMetaOnly:
    input: 
      "Combined_Analysis/super_contigs.{id}.fasta",
      expand("Combined_Analysis/subsample_superContig_alignmentReport.{subsample}.txt", subsample=subsampleIDs)
    output: "Combined_Analysis/super_contigs.{id}.alignment_report.txt"
    params:
      name="merge_superContig_alignment_miniMetaOnly",
      qos="normal",
      time="2:00:00",
      partition="normal",
      mem="4000", # don't change
      contig_thresh=parameters.ix['biosample_contig_thresh','entry']
    threads: 1
    version: "2.0"
    run: 
      # Managing files and obtain scratch location
      scratch = os.environ["LOCAL_SCRATCH"]
      input_on_scratch = names_on_scratch(input, scratch)
      output_on_scratch = names_on_scratch(output, scratch)
      contig_on_scratch = names_on_scratch(["super_contig_subset.fasta"], scratch)
      cp_to_scratch(input, scratch)
      # Perform organization of contigs    
      shell("""
        source activate {python2_env}
        python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[0]} {contig_on_scratch[0]}
        """)
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


rule make_supercontig_indices:
  input:
    "Combined_Analysis/super_contigs.{id}.fasta"
  output:
    "Combined_Analysis/superContigIndex_{id}.1.bt2l",
    "Combined_Analysis/superContigIndex_{id}.2.bt2l",
    "Combined_Analysis/superContigIndex_{id}.3.bt2l",
    "Combined_Analysis/superContigIndex_{id}.4.bt2l",
    "Combined_Analysis/superContigIndex_{id}.rev.1.bt2l",
    "Combined_Analysis/superContigIndex_{id}.rev.2.bt2l"
  params:
    name="make_supercontig_indices",
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
      basename=$(echo {output[0]} | cut -d/ -f2 | cut -d. -f1)
      echo $basename
      date; cd {scratch}; source activate {python2_env}
      python {code_dir}/process_scaffolds.py --lengthThresh {params.contig_thresh} {input_on_scratch} threshold_super_contigs.fasta
      source activate {python3_env}
      bowtie2-build --quiet --large-index --threads {threads} -f threshold_super_contigs.fasta $basename
      source deactivate; date
      """)
    # Move the rest of the reads back
    cp_from_scratch(output, scratch)



rule subsampleSuperContigPileup:
  input: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq",
    "{subsample}/S1.{subsample}.fastq",
    "{subsample}/S2.{subsample}.fastq",
    expand("Combined_Analysis/super_contigs.{id}.fasta", id=biosample),
    expand("Combined_Analysis/superContigIndex_{id}.1.bt2l", id=biosample),
    expand("Combined_Analysis/superContigIndex_{id}.2.bt2l", id=biosample),
    expand("Combined_Analysis/superContigIndex_{id}.3.bt2l", id=biosample),
    expand("Combined_Analysis/superContigIndex_{id}.4.bt2l", id=biosample),
    expand("Combined_Analysis/superContigIndex_{id}.rev.1.bt2l", id=biosample),
    expand("Combined_Analysis/superContigIndex_{id}.rev.2.bt2l", id=biosample)
    # expand("Combined_Analysis/super_contigs.{id}.fasta", id=biosample)
  output:
    temp("Combined_Analysis/subsample_superContig_mpileupReport.{subsample}.txt"),
    temp("Combined_Analysis/subsampleAlign2Supercontig.{subsample}.bam"),
    temp("Combined_Analysis/subsample_superContig_alignmentReport.{subsample}.txt")
  params:
    name="subsampleSuperContigPileup",
    qos="normal",
    time="23:59:59",
    partition=parameters.ix['subsample_bowtie2_partition','entry'],
    mem=parameters.ix['subsample_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Bowtie2 Alignment of reads back to super_contigs (contigs with new names)
    # This is perhaps not the best way, use the line above instead
    # awk 'a1==$1 {a2+=$4; next} {print a1,"\t", a2; a1=$1; a2=$4} END {print a1,"\t",a2}' alignResults.pile > {output_on_scratch[2]}
    shell("""
      basename=$(echo {input[5]} | cut -d/ -f2 | cut -d. -f1)
      echo {scratch}; cd {scratch}
      source activate {python2_env}
      python {code_dir}/process_scaffolds.py --lengthThresh {params.contig_thresh} {input_on_scratch[4]} threshold_super_contigs.fasta
      source activate {python3_env}
      cat {input_on_scratch[2]} {input_on_scratch[3]} > single.fastq
      bowtie2 --very-sensitive-local -I 100 -X 2000 -p {threads} -t -x $basename -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -U single.fastq -S alignResults.sam
      samtools view -b -o alignResults.bam alignResults.sam
      samtools sort -o alignResults_sorted.bam alignResults.bam
      cp alignResults_sorted.bam {output_on_scratch[1]}
      samtools index alignResults_sorted.bam
      samtools mpileup -f threshold_super_contigs.fasta -o alignResults.pile alignResults_sorted.bam 
      awk '$4>=5 {{print}}' alignResults.pile | cut -f 1 | sort | uniq -c | awk '{{print $2,"\t",$1}}' > {output_on_scratch[0]}
      # echo nextLineIsTheProblem
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
    qos="normal",
    time="1-0",
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
    contig_on_scratch = names_on_scratch(["super_contig_subset.fasta"], scratch)
    bamfiles_on_scratch = input_on_scratch[1:]
    cp_to_scratch(input, scratch)
    # Perform organization of contigs    
    shell("""
      echo {scratch}; cd {scratch} 
      source activate {python2_env}
      python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[0]} {contig_on_scratch[0]}
      source activate {python3_env}
      samtools mpileup -g -t DP -f {contig_on_scratch[0]} -o subsample_alignment.bcf {bamfiles_on_scratch}
      echo BCF_file_generation_completed
      samtools mpileup -f {contig_on_scratch[0]} -o {output_on_scratch[1]} {bamfiles_on_scratch}
      bcftools call -c -v -o {output_on_scratch[0]} subsample_alignment.bcf 
      echo VCF_file_generation_completed
      """)
    cp_from_scratch(output, scratch)
  

if bulk_flag =='Yes' or bulk_flag == 'yes' or bulk_flag == 'Y' or bulk_flag == 'y':
  rule subsampleGenomeSize_withBulk:
    input: 
      "Combined_Analysis/super_contigs.{id}.fasta",
      expand("Combined_Analysis/subsample_superContig_mpileupReport.{subsample}.txt", subsample=subsampleIDs),
      expand("Combined_Analysis/subsample_superContig_alignmentReport.{bulksample}.txt", bulksample=bulksampleIDs)
    output: 
      "Combined_Analysis/super_contigs.{id}.subsampleGenomeSize.txt"
    params:
      name="subsampleGenomeSize_withBulk",
      qos="normal",
      time="2:00:00",
      partition="normal",
      mem="4000", # don't change
      contig_thresh=parameters.ix['biosample_contig_thresh','entry']
    threads: 1
    version: "2.0"
    run: 
      # Managing files and obtain scratch location
      scratch = os.environ["LOCAL_SCRATCH"]
      input_on_scratch = names_on_scratch(input, scratch)
      output_on_scratch = names_on_scratch(output, scratch)
      contig_on_scratch = names_on_scratch(["super_contig_subset.fasta"], scratch)
      cp_to_scratch(input, scratch)
      # Perform organization of contigs    
      shell("""
        source activate {python2_env}
        python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[0]} {contig_on_scratch[0]}
        """)
      input_on_scratch[0] = contig_on_scratch[0]
      # What the next line does is to compute the NUMBER OF BPs COVERED (not coverage depth)
      merge_superContig_alignment_results(input_on_scratch, output_on_scratch)
      cp_from_scratch(output, scratch)
else:
  rule subsampleGenomeSize_miniMetaOnly:
    input: 
      "Combined_Analysis/super_contigs.{id}.fasta",
      expand("Combined_Analysis/subsample_superContig_mpileupReport.{subsample}.txt", subsample=subsampleIDs)
    output: 
      "Combined_Analysis/super_contigs.{id}.subsampleGenomeSize.txt"
    params:
      name="subsampleGenomeSize_miniMetaOnly",
      qos="normal",
      time="2:00:00",
      partition="normal",
      mem="4000", # don't change
      contig_thresh=parameters.ix['biosample_contig_thresh','entry']
    threads: 1
    version: "2.0"
    run: 
      # Managing files and obtain scratch location
      scratch = os.environ["LOCAL_SCRATCH"]
      input_on_scratch = names_on_scratch(input, scratch)
      output_on_scratch = names_on_scratch(output, scratch)
      contig_on_scratch = names_on_scratch(["super_contig_subset.fasta"], scratch)
      cp_to_scratch(input, scratch)
      # Perform organization of contigs    
      shell("""
        source activate {python2_env}
        python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[0]} {contig_on_scratch[0]}
        """)
      input_on_scratch[0] = contig_on_scratch[0]
      # What the next line does is to compute the NUMBER OF BPs COVERED (not coverage depth)
      merge_superContig_alignment_results(input_on_scratch, output_on_scratch)
      cp_from_scratch(output, scratch)

