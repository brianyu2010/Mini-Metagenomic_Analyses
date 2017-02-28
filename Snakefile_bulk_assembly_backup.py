###############################################
# Snakemake rules associated with assembly
# blast, quast of bulk subsamples.
# this file must be included into another 
# Snakemake file
#
# 2017.01.13 Added metaspade assembly in addition
#            megahit assembly
# 2017.02.01 Edited to work with Sherlock
###############################################

# rules

rule megahit_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2015.11.30 This rule is still using splitNcluster outputs.
  #            These need to be changed to the P1 P2 S1 S2 from trim
  # 2015.12.01 Changes were made to use trimmed reads directly
  # 2017.01.13 Updated with megahit functions for Bulk
  # 2017.02.03 Abundance files (not created) should be like this
  #            contig_1_name\tabundance
  #            contig_2_name\tabundance
  input: 
    "{subsample}/P1.{subsample}.fastq",
    "{subsample}/P2.{subsample}.fastq",
    "{subsample}/S1.{subsample}.fastq", 
    "{subsample}/S2.{subsample}.fastq"
  output: 
    "{subsample}/megahit_contigs.{subsample}.fasta",
    "{subsample}/quast_report_megahitBulk.{subsample}.txt"
  params: 
    name="bulk_megahit_assembly",
    qos="normal",
    time="2-0",
    partition=parameters.ix['bulk_megahit_partition','entry'], 
    mem=parameters.ix['bulk_megahit_memory','entry'],
    contig_thresh=parameters.ix['biosample_contig_thresh','entry'],
    kmer=parameters.ix['bulk_megahit_kmerlist','entry']
  threads: int(parameters.ix['bulk_megahit_thread','entry'])
  version: "3.0"
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
      # Assembly using megahit (for metagenomes). 
      shell("""
        # No longer cd to scratch
        echo {scratch}; date; pwd; echo
        source activate {python2_env}
        mem=$( echo {params.mem} )000000
        echo 'memory used in Bytes is '$mem 
        megahit_output_dir={wildcards.subsample}/megahit_output_{wildcards.subsample}
        source activate {python3_env}
        echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( cat {input[0]} | wc -l ) / 4 ))
        cat {input[2]} {input[3]} > {scratch}/single.fastq
        echo -e 'Starting Assembly: Number of single end reads is: \t'$(( $( wc -l < {scratch}/single.fastq) / 4 ))
        if [ -d $megahit_output_dir ]
        then
        time megahit --k-list {params.kmer} -t {threads} -m $mem --mem-flag 1 --continue -1 {input[0]} -2 {input[1]} -r {scratch}/single.fastq -o $megahit_output_dir --out-prefix megahit
        else
        time megahit --k-list {params.kmer} -t {threads} -m $mem --mem-flag 1 -1 {input[0]} -2 {input[1]} -r {scratch}/single.fastq -o $megahit_output_dir --out-prefix megahit
        fi
        date; ls {scratch}; echo; ls $megahit_output_dir
        cp $megahit_output_dir/megahit.contigs.fa {output_on_scratch[0]}
        echo 'Assembly Completed'
        source deactivate; pwd
        """)
      addToContigName(output_on_scratch[0], 0, wildcards.subsample, output[0])
      # cp_from_scratch(output[0], scratch)
      # Performing Quast
      shell("""
        echo 'Performing metaquast'
        date; pwd; echo
        source activate {python2_env}
        metaquast_output_dir={wildcards.subsample}/metaquast_megahit
        quast_output_dir={scratch}/quast_output
        metaquast.py --plots-format svg --gene-finding --max-ref-number 200 -m {params.contig_thresh} -t {threads} -o $metaquast_output_dir {output[0]}
        quast.py -m {params.contig_thresh} -t {threads} -o $quast_output_dir {output[0]}
        cp $quast_output_dir/report.txt {output[1]}
        source deactivate
        """)
      assert(file_empty(output)),"The contig file is empty."
    else:
      print('Input files have size 0')


rule metaSPAdes_assembly:
  # order of the inputs are paired1 paired2 single1 single2
  # 2017.01.13 Added this to compare with megahit assembly on bulk
  input: 
    "{subsample}/P1.{subsample}.fastq",
    "{subsample}/P2.{subsample}.fastq"
    # "{subsample}/S1.{subsample}.fastq", 
    # "{subsample}/S2.{subsample}.fastq"
  output: 
    "{subsample}/metaSPAdes_contigs.{subsample}.fasta",
    "{subsample}/quast_report_metaSPAdesBulk.{subsample}.txt"
  params: 
    name="bulk_metaSPAdes_assembly",
    qos="normal",
    time="3-0",
    partition=parameters.ix['biosample_assembly_partition','entry'], 
    mem=parameters.ix['biosample_assembly_memory','entry'],
    contig_thresh=parameters.ix['biosample_contig_thresh','entry'],
    kmer=parameters.ix['bulk_spades_kmerlist','entry']
  threads: int(parameters.ix['biosample_assembly_thread','entry'])
  version: "3.0"
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
      # Assembly using metaSPAdes (for shotgun metagenomes). 
      # now do assembly, keep output folder in scratch
      shell("""
        echo {scratch}; date; pwd; echo
        source activate {python2_env}
        mem=$( echo {params.mem} | rev | cut -c 4- | rev ) 
        echo 'memory used in Gb is '$mem 
        spades_output_dir={wildcards.subsample}/metaSPAdes_output_{wildcards.subsample}
        source activate {python3_env}
        echo -e 'Starting Assembly: Number of paired end reads is: \t'$(( $( cat {input[0]} | wc -l ) / 4 ))
        if [ -d $spades_output_dir ]
        then
        time metaspades.py -k {params.kmer} -t {threads} --continue -m $mem -1 {input[0]} -2 {input[1]} -o $spades_output_dir
        else
        time metaspades.py -k {params.kmer} -t {threads} -m $mem -1 {input[0]} -2 {input[1]} -o $spades_output_dir
        fi
        date; ls {scratch}; echo; ls $spades_output_dir
        cp $spades_output_dir/contigs.fasta {output_on_scratch[0]}
        echo 'Assembly Completed'
        source deactivate; pwd
        """)
      addToContigName(output_on_scratch[0], 0, wildcards.subsample, output[0])
      # cp_from_scratch(output[0], scratch)
      # Performing Quast
      shell("""
        echo 'Performing metaquast'
        date; pwd; echo
        source activate {python2_env}
        metaquast_output_dir={wildcards.subsample}/metaquast_metaSPAdes
        quast_output_dir={scratch}/quast_output
        metaquast.py --plots-format svg --gene-finding --max-ref-number 200 -m {params.contig_thresh} -t {threads} -o $metaquast_output_dir {output[0]}
        quast.py -m {params.contig_thresh} -t {threads} -o $quast_output_dir {output[0]}
        cp $quast_output_dir/report.txt {output[1]}
        source deactivate
        """)
      assert(file_empty(output)),"The contig file is empty."
    else:
      print('Input files have size 0')


rule align_to_bulk_megahit_assembly:
  # The output of this rule can be made more informative by including coverage of each contigs.
  input: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq",
    "{subsample}/S1.{subsample}.fastq",
    "{subsample}/S2.{subsample}.fastq",
    "{subsample}/megahit_contigs.{subsample}.fasta"
  output:
    "{subsample}/megahit_contig_alignment.{subsample}.pile"
    # "{subsample}/megahit_Abundance.{subsample}.cnt"
  params: 
    name="bulkMegahitAbundance",
    qos="normal",
    time="23:00:00",
    mem_per_core="4G",
    partition=parameters.ix['subsample_bowtie2_partition','entry'],
    mem=parameters.ix['subsample_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    # Bowtie2 Alignment of reads back to bulk megahit contigs
    shell("""
      echo {scratch}; date; echo
      source activate {python2_env}
      echo 'Thresholding contigs'
      python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input[4]} {scratch}/temp.fasta
      source activate {python3_env};
      echo 'Creating bowtie2 indices'
      bowtie2-build --quiet --threads {threads} -f {scratch}/temp.fasta {scratch}/bulkContigBase
      cat {input[2]} {input[3]} > {scratch}/single.fastq
      echo; ls {scratch}; echo; du -sh {scratch}; echo
      echo 'Aligning both paired and single reads'
      bowtie2 --phred33 --very-sensitive-local -I 100 -X 2000 -p {threads} -t -x {scratch}/bulkContigBase -1 {input[0]} -2 {input[1]} -U {scratch}/single.fastq -S {wildcards.subsample}/alignResults.sam
      echo; ls {scratch}; echo; du -sh {scratch}; echo;
      samtools_temp_dir={scratch}/temp_output/ # This is for samtools sort to put temp files
      samtools view -b -o {wildcards.subsample}/alignResults.bam {wildcards.subsample}/alignResults.sam
      samtools sort -m {params.mem_per_core} --threads {threads} -T $samtools_temp_dir -o {wildcards.subsample}/alignResults_sorted.bam {wildcards.subsample}/alignResults.bam
      samtools index {wildcards.subsample}/alignResults_sorted.bam
      samtools mpileup -f {scratch}/temp.fasta -o {output[0]} {wildcards.subsample}/alignResults_sorted.bam
      rm {wildcards.subsample}/alignResults* # Remove intermediate files in  directory
      echo 'Currently this rule only returns the pileup file'; date
      """)
    cp_from_scratch(output, scratch)


rule align_to_bulk_metaSPAdes_assembly:
  # The output of this rule can be made more informative by including coverage of each contigs.
  input: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq",
    "{subsample}/S1.{subsample}.fastq",
    "{subsample}/S2.{subsample}.fastq",
    "{subsample}/metaSPAdes_contigs.{subsample}.fasta"
  output:
    "{subsample}/metaSPAdes_contig_alignment.{subsample}.pile"
    # "{subsample}/metaSPAdes_Abundance.{subsample}.cnt"
  params: 
    name="bulkMetaSPAdesAbundance",
    qos="normal",
    time="23:00:00",
    mem_per_core="4G",
    partition=parameters.ix['subsample_bowtie2_partition','entry'],
    mem=parameters.ix['subsample_bowtie2_memory','entry'],
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: int(parameters.ix['subsample_bowtie2_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    # Bowtie2 Alignment of reads back to bulk megahit contigs
    shell("""
      echo {scratch}; date; echo
      source activate {python2_env}
      echo 'Thresholding contigs'
      python {code_dir}/threshold_scaffolds.py {params.contig_thresh} {input[4]} {scratch}/temp.fasta
      source activate {python3_env}
      echo 'Creating bowtie2 indices'
      bowtie2-build --quiet --threads {threads} -f {scratch}/temp.fasta {scratch}/bulkContigBase
      cat {input[2]} {input[3]} > {scratch}/single.fastq
      echo; ls {scratch}; echo; du -sh {scratch}; echo  
      echo 'Aligning both paired and single reads'
      bowtie2 --phred33 --very-sensitive-local -I 100 -X 2000 -p {threads} -t -x {scratch}/bulkContigBase -1 {input[0]} -2 {input[1]} -U {scratch}/single.fastq -S {wildcards.subsample}/alignResults.sam
      echo; ls {scratch}; echo; du -sh {scratch}; echo;
      samtools_temp_dir={scratch}/temp_output/ # This is for samtools sort to put temp files
      samtools view -b -o {wildcards.subsample}/alignResults.bam {wildcards.subsample}/alignResults.sam
      samtools sort -m {params.mem_per_core} --threads {threads} -T $samtools_temp_dir -o {wildcards.subsample}/alignResults_sorted.bam {wildcards.subsample}/alignResults.bam
      samtools index {wildcards.subsample}/alignResults_sorted.bam
      samtools mpileup -f {scratch}/temp.fasta -o {output[0]} {wildcards.subsample}/alignResults_sorted.bam
      echo 'Currently this rule only returns the pileup file'; date
      """)
    cp_from_scratch(output, scratch)

