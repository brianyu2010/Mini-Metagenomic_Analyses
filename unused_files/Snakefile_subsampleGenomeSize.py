###############################################
# Snakemake rules associated with analyzing 
# super contigs in conjunction with subsamples
# to compute size of genome from each subsample
# this file must be included into another 
# Snakemake file 
#
# 2016.08.30 Created Brian Yu
###############################################   

# In the future make this work better and then 
# comment out the rule subsampleReads_Align2SuperContigs
# in Snakefile_superContigAnalysis.py

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
    mem="20000", #parameters.ix['subsample_bowtie2_memory','entry'],
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
    shell("""cp {code_dir}/threshold_scaffolds.py {scratch}
      cd {scratch}; set +u
      source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
      python threshold_scaffolds.py {params.contig_thresh} {input_on_scratch[2]} temp.fasta
      {tool_dir}/bowtie2-2.2.6/bowtie2-build -f temp.fasta spadeContigs
      {tool_dir}/bowtie2-2.2.6/bowtie2 --very-sensitive-local -I 0 -X 1000 -p {threads} -t -x spadeContigs -1 {input_on_scratch[0]} -2 {input_on_scratch[1]} -S alignResults.sam
      {tool_dir}/samtools-1.3/samtools view -b -o alignResults.bam alignResults.sam
      {tool_dir}/samtools-1.3/samtools sort -o alignResults_sorted.bam alignResults.bam
      cp alignResults_sorted.bam {output_on_scratch[1]}
      {tool_dir}/samtools-1.3/samtools index alignResults_sorted.bam
      {tool_dir}/samtools-1.3/samtools mpileup -f temp.fasta -o alignResults.pile alignResults_sorted.bam 
      awk '$4>=5 {{print}}' alignResults.pile | cut -f 1 | sort | uniq -c | awk '{{print $2,"\t",$1}}' > {output_on_scratch[0]}
      awk 'a1==$1 {a2+=$4; next} {print a1,"\t", a2; a1=$1; a2=$4} END {print a1,"\t",a2}' alignResults.pile > {output_on_scratch[2]}
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
    mem="50000",
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: 1
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

