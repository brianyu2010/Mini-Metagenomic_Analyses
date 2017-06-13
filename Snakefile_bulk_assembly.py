###############################################
# Snakemake rules associated with assembly
# blast, quast of bulk subsamples.
# this file must be included into another 
# Snakemake file
#
# 2017.01.13 Added metaspade assembly in addition
#            megahit assembly
# 2017.02.01 Edited to work with Sherlock
# 2017.02.24 Added bulk concatenate and quality trim rules
###############################################

# rules

# use an input file to specify and cat the fastq files into local10G
rule BulkConcatenate:
  # for each subsample, find fastq files from all sequencing runs and combine them.
  # need to be able to handle certain subsamples appearing on in some sequencing runs but not others
  # input: No Input
  # temp removes the files once it has been used
  output:
    temp("{subsample}/read1_bulk.{subsample}.fastq"),
    temp("{subsample}/read2_bulk.{subsample}.fastq")
  params: 
    name="bulk_combine_fastq",
    qos="normal",
    time="12:00:00",
    partition="normal,quake", 
    mem="125000" # Don't change this
  threads: 10
  version: "2.0"
  run: 
    scratch = os.environ["LOCAL_SCRATCH"]
    output_on_scratch = names_on_scratch(output, scratch)
    print(scratch)
    # 2017.02.02 Added to solve the ambiguity of subsample and bulksample
    # 2017.02.24 Changed to look at only bulk subsamples
    if bulk_flag:
      s = bulk_table.ix[wildcards.subsample,:]
      # if only one row, automatically become a column so need to handle that
      if len(s.shape) == 1:
        print('Sample '+wildcards.subsample+' exists in only one sequencing run.')
        fastq_dir = s.ix["biosample"] + "/" + s.ix["sequencingrun"] + "/" + s.ix["subsamplename"]
        shell("""
          cd {fastq_dir}; ls
          cp *.fastq.gz {scratch}
          cd {scratch}; echo; du -sh; echo
            ls *.fastq.gz > zipped_file_names.txt
          while read line; do gzip -d $line; done < zipped_file_names.txt
          echo; du -sh; echo
          ls *_R1_001.fastq | sort > read1files.txt; head read1files.txt
          ls *_R2_001.fastq | sort > read2files.txt; head read2files.txt
          """)
        # Putting in a separate shell script so that pwd is not on scratch
        shell("""
          echo; pwd; date; echo
          while read line; do cat {scratch}/$line >> {output[0]}
          echo -e 'Number of reads in read1.fastq is: '$(( $( wc -l < {scratch}/$line ) / 4 ))
          rm {scratch}/$line; done < {scratch}/read1files.txt
          echo 'Finished read 1'
          while read line; do cat {scratch}/$line >> {output[1]}
          echo -e 'Number of reads in read2.fastq is: '$(( $( wc -l < {scratch}/$line ) / 4 ))
          rm {scratch}/$line; done < {scratch}/read2files.txt
          echo
          echo -e 'Number of reads in Read1.output is: '$(( $( wc -l < {output[0]}) / 4 ))
          echo -e 'Number of reads in Read2.output is: '$(( $( wc -l < {output[1]}) / 4 ))
          echo
          echo 'Combining fastq files completed'; date
          """)
      else:
        numseqruns = len(s.index)
        print('Sample '+wildcards.subsample+' exists in '+str(numseqruns)+' sequencing runs.')
        for i in range(numseqruns):
          fastq_dir = s.ix[i,"biosample"] + "/" + s.ix[i,"sequencingrun"] + "/" + s.ix[i,"subsamplename"]
          shell("""
            cd {fastq_dir}; ls
            cp *.fastq.gz {scratch}
            cd {scratch}; echo; du -sh; echo
            ls *.fastq.gz > zipped_file_names.txt
            while read line; do gzip -d $line; done < zipped_file_names.txt
            echo; du -sh; echo
            ls *_R1_001.fastq | sort > read1files.txt; head read1files.txt
            ls *_R2_001.fastq | sort > read2files.txt; head read2files.txt
            """)
          # Putting in a separate shell script so that pwd is not on scratch
          shell("""
            echo; pwd; date; echo
            while read line; do cat {scratch}/$line >> {output[0]}
            echo -e 'Number of reads in read1.fastq is: '$(( $( wc -l < {scratch}/$line ) / 4 ))
            rm {scratch}/$line; done < {scratch}/read1files.txt
            while read line; do cat {scratch}/$line >> {output[1]}
            echo -e 'Number of reads in read2.fastq is: '$(( $( wc -l < {scratch}/$line ) / 4 ))
            rm {scratch}/$line; done < {scratch}/read2files.txt
            echo
            echo -e 'Number of reads in Read1.output is: '$(( $( wc -l < {output[0]}) / 4 ))
            echo -e 'Number of reads in Read2.output is: '$(( $( wc -l < {output[1]}) / 4 ))
            echo
            echo 'Combining fastq files completed'; date
            """)
      assert(file_empty(output)),"Output fastq files are empty."
      assert(check_lines(output[0],output[1])),"Output fastq files have different number of lines."
      assert(check_fastq_ids(output[0],output[1])),"Output fastq ids are not in the same order."
    else:
      print('Error: bulk sample flag is not asserted')

    

rule bulk_quality_trim:
  input: 
    "{subsample}/read1_bulk.{subsample}.fastq", 
    "{subsample}/read2_bulk.{subsample}.fastq"
  output: 
    "{subsample}/P1_bulk.{subsample}.fastq", 
    "{subsample}/S1_bulk.{subsample}.fastq", 
    "{subsample}/P2_bulk.{subsample}.fastq", 
    "{subsample}/S2_bulk.{subsample}.fastq"
  params: 
    name="bulk_quality_trim", 
    qos="normal",
    time="12:00:00",
    partition="normal,quake", 
    mem="125000",
    trim_to_read_length=str(parameters.ix['Desired_Read_Length','entry']),
    downsample_read_number=str(4*int(parameters.ix['Down_Sample_Read_Number','entry']))
  threads: 10
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["LOCAL_SCRATCH"]
    input_on_scratch = names_on_scratch(input, scratch)
    # output_on_scratch = names_on_scratch(output, scratch)
    # cp_to_scratch(input, scratch) # No need
    # trim fastq to desired read length (ie. 75 bp) if required
    if int(params.trim_to_read_length) <= 0:
      print('Do no trim reads')
      shell("""cp {input[0]} {scratch}/trimmed1.fastq; cp {input[1]} {scratch}/trimmed2.fastq""")
    else:
      shell("""
        source activate {python3_env}
        fastx_trimmer -Q33 -l {params.trim_to_read_length} -i {input[0]} -o {scratch}/trimmed1.fastq}
        fastx_trimmer -Q33 -l {params.trim_to_read_length} -i {input[1]} -o {scratch}/trimmed2.fastq}
        """)
    # randomize fastq using a python script
    if int(params.downsample_read_number) <= 0:
      print('Do not downsample reads')
      shell("""mv {scratch}/trimmed1.fastq {input_on_scratch[0]}
        mv {scratch}/trimmed2.fastq {input_on_scratch[1]}""")
    else:
      # take the frist set number of lines
      print('Down Sample Read Number of Lines is: '+params.downsample_read_number)
      shell("""
        source activate {python2_env}
        python --version
        python {code_dir}/randomizeFastq.py {scratch}/trimmed1.fastq {scratch}/trimmed2.fastq {scratch}/temp
        head -n {params.downsample_read_number} {scratch}/temp1.fastq > {input_on_scratch[0]}
        head -n {params.downsample_read_number} {scratch}/temp2.fastq > {input_on_scratch[1]}
        wc -l {input_on_scratch[0]}
        wc -l {input_on_scratch[1]}
        rm {scratch}/temp*.fastq
        """)
    # Performing trimmomatic trimming; note there are 4 outputs
    # Trimmomatic can be multi-threaded
    shell("""
      echo {scratch}; date; echo; du -sh {scratch}; echo
      source activate {python3_env}
      cp {code_dir}/Combined_PE_V2.fa {scratch}/adapterSeqs.fa
      default_seedMismatch=3   # not sensitive
      default_palen_th=30      # not sensitive
      default_minAdapterLen=3  # not sensitive
      default_slidingWindow=10 # this is a sensitive parameter. change to 6 and you get 50% less reads
      default_slidingQual=30   # I want this sliding window quality of 30
      default_maxInfoLen=120   # This is your target read length, possible choices are 120, 125, 140
      default_maxInfo_th=0.5   # Not very sensitive to this either. Probably because Q30 is dominant
      default_option_string="ILLUMINACLIP:{scratch}/adapterSeqs.fa:$default_seedMismatch:$default_palen_th:10:$default_minAdapterLen:TRUE SLIDINGWINDOW:$default_slidingWindow:$default_slidingQual MAXINFO:$default_maxInfoLen:$default_maxInfo_th LEADING:30 TRAILING:30 MINLEN:30"
      echo 'Start Trimming with Trimmomatic'
      trimmomatic PE -phred33 -threads {threads} -trimlog {scratch}/pairtrim.log {input_on_scratch} {output} $default_option_string
      echo 'Trimming Completed'; date; source deactivate
      """)
    assert(check_lines(output[0],output[2])),"Paired output files have different number of lines."
    assert(check_fastq_ids(output[0],output[2])),"Paired output files have different read id numbers/orders."


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
    "{subsample}/P1_bulk.{subsample}.fastq",
    "{subsample}/P2_bulk.{subsample}.fastq",
    "{subsample}/S1_bulk.{subsample}.fastq", 
    "{subsample}/S2_bulk.{subsample}.fastq"
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
    "{subsample}/P1_bulk.{subsample}.fastq",
    "{subsample}/P2_bulk.{subsample}.fastq"
    # "{subsample}/S1_bulk.{subsample}.fastq", 
    # "{subsample}/S2_bulk.{subsample}.fastq"
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

      
"""
rule align_to_bulk_megahit_assembly:
  # The output of this rule can be made more informative by including coverage of each contigs.
  input: 
    "{subsample}/P1_bulk.{subsample}.fastq", 
    "{subsample}/P2_bulk.{subsample}.fastq",
    "{subsample}/S1_bulk.{subsample}.fastq",
    "{subsample}/S2_bulk.{subsample}.fastq",
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
    shell(""
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
      "")
    cp_from_scratch(output, scratch)


rule align_to_bulk_metaSPAdes_assembly:
  # The output of this rule can be made more informative by including coverage of each contigs.
  input: 
    "{subsample}/P1_bulk.{subsample}.fastq", 
    "{subsample}/P2_bulk.{subsample}.fastq",
    "{subsample}/S1_bulk.{subsample}.fastq",
    "{subsample}/S2_bulk.{subsample}.fastq",
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
    shell(""
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
      "")
    cp_from_scratch(output, scratch)
"""
