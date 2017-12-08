#######################################################
# Rules to align all subsamples and shotgun samples
# back to the re-assembled genomes.
#
# Added 2017.10.20
#######################################################

rule make_reassembled_grouped_genome_index:
  input:
    "Genome_Reassembly/genome_scaffolds_withBulk.{genome}.fasta"
  output:
    temp("Genome_Reassembly/genome_scaffolds_withBulk_{genome}.1.bt2l"),
    temp("Genome_Reassembly/genome_scaffolds_withBulk_{genome}.2.bt2l"),
    temp("Genome_Reassembly/genome_scaffolds_withBulk_{genome}.3.bt2l"),
    temp("Genome_Reassembly/genome_scaffolds_withBulk_{genome}.4.bt2l"),
    temp("Genome_Reassembly/genome_scaffolds_withBulk_{genome}.rev.1.bt2l"),
    temp("Genome_Reassembly/genome_scaffolds_withBulk_{genome}.rev.2.bt2l")
  params:
    name="reassembled_grouped_genome_index",
    qos="normal",
    time="12:00:00",
    partition="quake,normal",
    mem="80000",
    contig_thresh=parameters.ix['reassembly_contig_thresh','entry']
  threads: 5
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing bowtie2-build
    shell("""
      basename=$(echo {output[0]} | cut -d/ -f2 | cut -d. -f1)
      echo $basename
      date; cd {scratch};
      # Threshold scaffolds
      source activate {python2_env}
      python {code_dir}/process_scaffolds.py --lengthThresh {params.contig_thresh} {input_on_scratch} threshold_super_contigs.fasta
      source activate {python3_env}
      bowtie2-build --quiet --large-index --threads {threads} -f threshold_super_contigs.fasta $basename
      source deactivate; date
      """)
    # Move the rest of the reads back
    cp_from_scratch(output, scratch)


rule reassembledGenome_minimeta_alignment_allSubSamples:
  input:
    expand("{subsample}/P1.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/P2.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/S1.{subsample}.fastq", subsample=subsampleIDs),
    expand("{subsample}/S2.{subsample}.fastq", subsample=subsampleIDs),
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.1.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.2.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.3.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.4.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.rev.1.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.rev.2.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk.{genome}.fasta"
  output:
    # expand("{subsample}/alignResults_{{genome}}.bam",subsample=subsampleIDs) 
    "Genome_Reassembly/miniMetaReads_realignmentDepthProfile.{genome}.txt"
  params:
    name="reassembledGenome_minimeta_alignment",
    qos="normal",
    time=parameters.ix['reassembly_bowtie2_time','entry'],
    mem_per_core="10G",
    partition=parameters.ix['reassembly_bowtie2_partition','entry'],
    mem=parameters.ix['reassembly_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['reassembly_contig_thresh','entry']
  threads: int(parameters.ix['reassembly_bowtie2_thread','entry'])
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["L_SCRATCH_JOB"]
    output_on_scratch = names_on_scratch(output, scratch)
    # Performing bowtie2 alignment --time --phred33 --un <path> --un-conc <path>
    numSubSamples = len(subsampleIDs)
    p1 = input[(0 * numSubSamples) : (1 * numSubSamples)]
    p2 = input[(1 * numSubSamples) : (2 * numSubSamples)]
    s1 = input[(2 * numSubSamples) : (3 * numSubSamples)]
    s2 = input[(3 * numSubSamples) : (4 * numSubSamples)]
    # For each subsamples, create a bam file
    basename = input[-2].split('.')[0] # Last file is grouped scaffolds
    print(basename)
    for i in range(numSubSamples):
      sample = subsampleIDs[i]
      shell("echo {sample}")
      p1fastq = p1[i]
      p2fastq = p2[i]
      s1fastq = s1[i]
      s2fastq = s2[i]
      shell("""
        # no longer cd to scratch
        source activate {python3_env}
        cat {s1fastq} {s2fastq} > {scratch}/single_reads.{sample}.fastq
        date; echo
        bowtie2 --phred33 --very-sensitive -I 100 -X 2000 -p {threads} -x {basename} -1 {p1fastq} -2 {p2fastq} -U {scratch}/single_reads.{sample}.fastq -S {sample}/alignResults_{wildcards.genome}.sam
        echo; du -sh {scratch}; echo
        samtools_temp_dir={scratch}/temp_output/ # This is for samtools sort to put temp files
        samtools view -f 0x003 -bh -o {sample}/alignResults_{wildcards.genome}.bam {sample}/alignResults_{wildcards.genome}.sam
        samtools sort -m {params.mem_per_core} --threads {threads} -o {sample}/alignResults_{wildcards.genome}_sorted.bam {sample}/alignResults_{wildcards.genome}.bam
        samtools index {sample}/alignResults_{wildcards.genome}_sorted.bam
        source deactivate
        rm {scratch}/single_reads.{sample}.fastq
        """)
    # Combine coverage profile from all bams
    genome_sorted_bam_files = [x+'/alignResults_'+wildcards.genome+'_sorted.bam' for x in subsampleIDs]
    header_text = 'ContigName Position ' + ' '.join(subsampleIDs) # separated by space, later changed to \t 
    # header_text = 'ContigName\tPosition\t' + '\t'.join(subsampleIDs)
    print(header_text)
    shell("""
      echo Combining_coverage_into_depth_file
      echo {header_text} | tr ' ' '\t' > {output_on_scratch} 
      source activate {python3_env}
      ls {scratch}
      samtools depth -a {genome_sorted_bam_files} >> {output_on_scratch}
      head {output_on_scratch} 
      source deactivate
      """)
    # Remove intermediate files saved in each subsample folder
    for i in range(numSubSamples):
      sample = subsampleIDs[i]
      shell("rm {sample}/alignResults_{wildcards.genome}* ")
    cp_from_scratch(output, scratch)


rule combine_miniMetaRealignment_to_Reassembly:
  input:
    expand("{subsample}/alignResults_{{genome}}.bam",subsample=subsampleIDs) 
  output:
    # "Genome_Reassembly/miniMetaReads_realignmentDepthProfile.{genome}.txt"
  params:
    name="minimeat_combine_realignment",
    qos="normal",
    time=parameters.ix['reassembly_bowtie2_time','entry'],
    mem_per_core="10G",
    partition="bigmem",  # parameters.ix['reassembly_bowtie2_partition','entry'],
    mem="500000" # parameters.ix['reassembly_bowtie2_memory','entry'],
    # contig_thresh=parameters.ix['reassembly_contig_thresh','entry']
  threads: 2 # int(parameters.ix['reassembly_bowtie2_thread','entry'])
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["L_SCRATCH_JOB"]
    output_on_scratch = names_on_scratch(output, scratch)
    # Performing bowtie2 alignment --time --phred33 --un <path> --un-conc <path>
    numSubSamples = len(subsampleIDs)
    # Combine coverage profile from all bams
    header_text = 'ContigName Position ' + ' '.join(subsampleIDs) # separated by space, later changed to \t 
    print(header_text)
    shell("""
      echo Combining_coverage_into_depth_file
      echo {header_text} | tr ' ' '\t' > {output_on_scratch} 
      source activate {python3_env}
      ls {scratch}
      samtools depth -a {input} >> {output_on_scratch}
      head {output_on_scratch} 
      source deactivate
      """)
    # Remove intermediate files saved in each subsample folder
    for i in range(numSubSamples):
      sample = subsampleIDs[i]
      shell("rm {sample}/alignResults_{wildcards.genome}* ")
    cp_from_scratch(output, scratch)


rule reassembledGenome_bulk_alignment:
  input:
    expand("{bulksample}/P1.{bulksample}.fastq", bulksample=bulksampleIDs),
    expand("{bulksample}/P2.{bulksample}.fastq", bulksample=bulksampleIDs),
    expand("{bulksample}/S1.{bulksample}.fastq", bulksample=bulksampleIDs),
    expand("{bulksample}/S2.{bulksample}.fastq", bulksample=bulksampleIDs),
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.1.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.2.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.3.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.4.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.rev.1.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk_{genome}.rev.2.bt2l",
    "Genome_Reassembly/genome_scaffolds_withBulk.{genome}.fasta"
  output: 
    "Genome_Reassembly/shotgunReads_realignmentDepthProfile.{genome}.txt",
    "Genome_Reassembly/shotgunReads_absoluteRealignmentReadCount.{genome}.txt",
    "Genome_Reassembly/shotgunReads_normalizedRealignmentReadCount.{genome}.txt"
  params:
    name="reassembledGenome_bulk_alignment",
    qos="normal",
    time=parameters.ix['reassembly_bowtie2_time','entry'],
    mem_per_core="10G",
    partition=parameters.ix['reassembly_bowtie2_partition','entry'],
    mem=parameters.ix['reassembly_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['reassembly_contig_thresh','entry']
  threads: int(parameters.ix['reassembly_bowtie2_thread','entry'])
  version: "1.0"
  run:
    # Manage files and obtain scratch location
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    # Performing bowtie2 alignment --time --phred33 --un <path> --un-conc <path>
    numBulkSamples = len(bulksampleIDs)
    p1 = input[(0 * numBulkSamples) : (1 * numBulkSamples)]
    p2 = input[(1 * numBulkSamples) : (2 * numBulkSamples)]
    s1 = input[(2 * numBulkSamples) : (3 * numBulkSamples)]
    s2 = input[(3 * numBulkSamples) : (4 * numBulkSamples)]
    # For each shotgun sample, create a bam file
    basename = input[-2].split('.')[0] # Last input is grouped scaffolds
    print(basename)
    for i in range(numBulkSamples):
      sample = bulksampleIDs[i]
      shell("echo {sample}")
      p1fastq = p1[i]
      p2fastq = p2[i]
      s1fastq = s1[i]
      s2fastq = s2[i]
      shell("""
        # no longer cd to scratch
        source activate {python3_env}
        cat {s1fastq} {s2fastq} > {scratch}/single_reads.{sample}.fastq
        date; echo
        bowtie2 --phred33 --very-sensitive -I 100 -X 2000 -p {threads} -x {basename} -1 {p1fastq} -2 {p2fastq} -U {scratch}/single_reads.{sample}.fastq -S {sample}/alignResults_{wildcards.genome}.sam
        echo; du -sh {scratch}
        samtools_temp_dir={scratch}/temp_output/ # This is for samtools sort to put temp files
        samtools view -f 0x003 -bh -o {sample}/alignResults_{wildcards.genome}.bam {sample}/alignResults_{wildcards.genome}.sam
        samtools sort -m {params.mem_per_core} --threads {threads} -o {sample}/alignResults_{wildcards.genome}_sorted.bam {sample}/alignResults_{wildcards.genome}.bam
        samtools index {sample}/alignResults_{wildcards.genome}_sorted.bam 
        samtools idxstats {sample}/alignResults_{wildcards.genome}_sorted.bam | cut -f 1,3 | head -n -1 > {sample}/alignResults_{wildcards.genome}_coverage.{sample}.txt 
        source deactivate
        echo
        """)
    # Combine coverage profile from all bulk bams
    genome_sorted_bam_files = [x+'/alignResults_'+wildcards.genome+'_sorted.bam' for x in bulksampleIDs]
    header_text = 'ContigName Position ' + ' '.join(bulksampleIDs) # separated by space, later changed to \t
    print(header_text)
    shell("""
      echo Combining_coverage_into_depth_file
      echo {header_text} | tr ' ' '\t' > {output_on_scratch[0]}
      source activate {python3_env}
      samtools depth -a {genome_sorted_bam_files} >> {output_on_scratch[0]}
      head {output_on_scratch[0]}
      source deactivate
      """)
    # Threashold the scaffold file to combine alignment read counts
    scaffold_file = input[-1]
    scaffold_on_scratch = input_on_scratch[-1]
    shell("""
      source activate {python2_env}
      python {code_dir}/process_scaffolds.py --lengthThresh {params.contig_thresh} {scaffold_file} {scaffold_on_scratch}
      source deactivate
      """)
    # First element in input needs to be the contig file, which is the last element of variable input
    merge_superContig_alignment_results([input_on_scratch[-1]] + [x+'/alignResults_'+wildcards.genome+'_coverage.'+x+'.txt' for x in bulksampleIDs], [output_on_scratch[1]], False)
    merge_superContig_alignment_results([input_on_scratch[-1]] + [x+'/alignResults_'+wildcards.genome+'_coverage.'+x+'.txt' for x in bulksampleIDs], [output_on_scratch[2]], True)
    # Remove intermediate files saved in each shotgun sample folder
    for i in range(numBulkSamples):
      sample = bulksampleIDs[i]
      shell("rm {sample}/alignResults_{wildcards.genome}* ")
    cp_from_scratch(output, scratch)


def merge_superContig_alignment_results(input, output, normalize_flag):
  """
  input: a list of files in the form of subsample_superContig_alignmentReport.ILxxxxxx.txt
  output: table containing total coverage average coverage of each contig in each sample
          normalized by contig length
  normalize_flag: flag to normalize coverage (length of contig covered) by contig length
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
          if normalize_flag:
            superContig_coverage[l[0]][ind] = float(l[1])/float(superContig_length[l[0]])
          else:
            superContig_coverage[l[0]][ind] = float(l[1]) # 2017.03.31 removed length normalization
  # Processing output file. output should be a list of only one string
  # print(subsamples)
  with open(output[0],'w') as f:
    subsamples.insert(0, '')
    t = f.write('\t'.join(subsamples) + '\n')
    for key in superContig_coverage:
      # add the supercontig name in front of all the coverages
      # print(key + '\t' + '\t'.join([str(x) for x in superContig_coverage[key]]))
      # separate by tab and add new line at the end (need to turn int list to str list first)
      t = f.write(key + '\t' + '\t'.join([str(x) for x in superContig_coverage[key]]) + '\n')

