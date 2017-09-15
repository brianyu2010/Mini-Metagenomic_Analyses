##########################################################
# Rules for importing, trimming, and clustering fastq
##########################################################

################################################
# rules
#
# 2016.09.21 Added prepocess step to cut 
#            reads to certain length and 
#            downsample reads if desired.
# 2017.01.31 Integrated some of the processes
#            into the master script and removed
#            unnecessary .sh and .py files. Also
#            adapted for sherlock
#################################################

# use an input file to specify and cat the fastq files into local10G
rule concatenate:
  # for each subsample, find fastq files from all sequencing runs and combine them.
  # need to be able to handle certain subsamples appearing on in some sequencing runs but not others
  # input: No Input
  # temp removes the files once it has been used
  output:
    temp("{subsample}/read1.{subsample}.fastq"),
    temp("{subsample}/read2.{subsample}.fastq")
  params: 
    name="combine_fastq",
    qos="normal",
    time="6:00:00",
    partition="normal,quake,owners", 
    mem="25000", # Don't change this
  threads: 2
  version: "2.0"
  run: 
    scratch = os.environ["L_SCRATCH_JOB"]
    output_on_scratch = names_on_scratch(output, scratch)
    print(scratch)
    # 2017.02.02 Added to solve the ambiguity of subsample and bulksample
    # 2017.02.24 Edited to only handle the mini-metagenomic samples
    """
    if bulk_flag:
      if wildcards.subsample in sample_table.index.tolist():
        s = sample_table.ix[wildcards.subsample,:]
      else:
        s = bulk_table.ix[wildcards.subsample,:]
    else:
    """
    s = sample_table.ix[wildcards.subsample,:] 
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


rule quality_trim:
  input: 
    "{subsample}/read1.{subsample}.fastq", 
    "{subsample}/read2.{subsample}.fastq"
  output: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/S1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq", 
    "{subsample}/S2.{subsample}.fastq"
  params: 
    name="quality_trim", 
    qos="normal",
    time="12:00:00",
    partition=parameters.ix['subsample_trim_partition','entry'], 
    mem=parameters.ix['subsample_trim_memory','entry'],
    trim_to_read_length=str(parameters.ix['Desired_Read_Length','entry']),
    downsample_read_number=str(4*int(parameters.ix['Down_Sample_Read_Number','entry']))
  threads: int(parameters.ix['subsample_trim_thread','entry'])
  version: "2.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["L_SCRATCH_JOB"]
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
    # cp_from_scratch(output, scratch)


rule overrep_seq:
  # use "wild card" to specify for multiple parallel processes
  input: "{subsample}/{type}.{subsample}.fastq"
  output: "{subsample}/{type}.{subsample}.fastqc_results.txt"
  params: 
    name="fastqc",
    qos="normal",
    time="30:00",
    partition="normal,quake,owners", 
    mem="8000" # Don't change this
  threads: 2
  version: "2.0"
  run: 
    # Managing files and obtain scratch location
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Quality checks
    # java is assumed in the path, which comes with the anaconda3/ environment
    shell("""
      echo {scratch}; cd {scratch}; date
      source activate {python3_env}
      foldername=$( echo {input_on_scratch} | rev | sed s/./_/6 | rev )"c"
      echo $foldername
      fastqc --extract -t {threads} {input_on_scratch}
      mv $foldername/fastqc_data.txt {output_on_scratch}
      """)
    cp_from_scratch(output, scratch)

    
"""
# could combined this rule into the trimming rule
rule preprocess:
  input: 
    "{subsample}/read1.{subsample}.fastq",
    "{subsample}/read2.{subsample}.fastq"
  output:
    temp("{subsample}/read1_preprocess.{subsample}.fastq"),
    temp("{subsample}/read2_preprocess.{subsample}.fastq")
  params:
    name="reads_preprocess",
    qos="normal",
    time="2:00:00",
    partition="normal",
    mem="12000",
    trim_to_read_length=str(parameters.ix['Desired_Read_Length','entry']),
    downsample_read_number=str(4*int(parameters.ix['Down_Sample_Read_Number','entry']))
  threads: 3
  version: "2.0"
  run:
    # organize input
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    # trim fastq to desired read length (ie. 75 bp) if required
    if int(params.trim_to_read_length) <= 0:
        print('Do no trim reads')
        shell("cp {input[0]} {input_on_scratch[0]}; cp {input[1]} {input_on_scratch[1]}")
    else:
      shell(""
        source activate {python3_env}
        fastx_trimmer -Q33 -l {params.trim_to_read_length} -i {input[0]} -o {input_on_scratch[0]}
        fastx_trimmer -Q33 -l {params.trim_to_read_length} -i {input[1]} -o {input_on_scratch[1]}
        "")
    # randomize fastq using a python script
    if int(params.downsample_read_number) <= 0:
      print('Do not downsample reads')
      shell("mv {input_on_scratch[0]} {output[0]}; mv {input_on_scratch[1]} {output[1]}")
    else:
      # take the frist set number of lines
      print('Down Sample Read Number of Lines is: '+params.downsample_read_number)
      shell(""
        source activate {python2_env}
        python --version
        python {code_dir}/randomizeFastq.py {input_on_scratch} {scratch}/temp
        head -n {params.downsample_read_number} {scratch}/temp1.fastq > {output[0]}
        head -n {params.downsample_read_number} {scratch}/temp2.fastq > {output[1]}
        wc -l {output[0]}
        wc -l {output[1]}
        "")
"""

