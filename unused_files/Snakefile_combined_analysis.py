###############################################
# Snakemake rules associated with combined
# analysis of each biosample from subsample resutls
# this file must be included into another 
# Snakemake file
###############################################
"""
# rename each contig with subsample information
# this rule is no longer used 2015.07.15
rule combine_contigs:
  input: expand("{subsample}/contigs.{subsample}.fasta", subsample=subsampleIDs)
  output: "Combined_Analysis/subsample_contigs.{id}.fasta", "Combined_Analysis/subsample_contigs_name.{id}.txt"
  params:    
    name="combine_contigs", 
    partition="general", 
    mem="3000"
  threads: 1
  version: "1.0"
  run: 
    assert(file_empty(input)),"One of the input contigs is empty."
    shell("source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/ &&\
      python {code_dir}/snakehelper_combine_subsample_contigs.py {input} -o {output[0]} -l {output[1]}")
    assert(file_empty(output)),"Either the combined contigs or names is empty."


# this rule is no longer used 2015.07.15
rule superContig_distribution:
  # it is important that the order for input is query and then contig
  input: 
    "{folder}/super_contigs.{id}.fasta", 
    "{folder}/subsample_contigs.{id}.fasta",
    "{folder}/super_contigs_name.{id}.txt", 
    "{folder}/subsample_contigs_name.{id}.txt"
  output: 
    "{folder}/super_contigs_blast_report.{id}.txt", 
    "{folder}/super_contigs_distribution.{id}.txt"
  params:    
    name="superContig_distribution", 
    partition="general", 
    mem="20000", # don't change this
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: 4
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing Blast and Data Summary
    shell("bash {code_dir}/snakehelper_localblast.sh {scratch} {input_on_scratch[0]} {input_on_scratch[1]} {output_on_scratch[0]} {tool_dir} {code_dir} {threads} {params.contig_thresh}")
    # path on scratch contains absolute path information
    shell("source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/ &&\
      python {code_dir}/snakehelper_contig_similarity.py {input_on_scratch[2]} {input_on_scratch[3]} {output_on_scratch[0]} {output_on_scratch[1]}")
    cp_from_scratch(output, scratch)


# This rule is no longer used 2015.07.15
rule contig_similarity:
  input: 
    "{folder}/{filename}_contigs.{id}.fasta", 
    "{folder}/{filename}_contigs_name.{id}.txt"
  output: 
    "{folder}/{filename}_contigs.{id}.similarity_blast_report.txt", 
    "{folder}/{filename}_contigs.{id}.similarity_matrix.txt"
  params:    
    name="contig_similarity", 
    partition="long", 
    mem="40000",
    contig_thresh=parameters.ix['biosample_contig_thresh','entry']
  threads: 11 
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing Blast and Data Summary
    shell("bash {code_dir}/snakehelper_localblast.sh {scratch} {input_on_scratch[0]} {input_on_scratch[0]} {output_on_scratch[0]} {tool_dir} {code_dir} {threads} {params.contig_thresh}") # typically use {contig_thresh}
    shell("source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/ &&\
      python {code_dir}/snakehelper_contig_similarity.py {input_on_scratch[1]} {input_on_scratch[1]} {output_on_scratch[0]} {output_on_scratch[1]}")
    assert(file_empty(output_on_scratch)),"One of the output files are empty."
    cp_from_scratch(output, scratch)



rule organize_subsample_blast:
  input: expand("{subsample}/BlastResults.{subsample}.txt", subsample=subsampleIDs)
  output: "Combined_Analysis/subsample_species_abundance.{id}.txt"
  params:    
    name="organize_subsample_blast", 
    partition="general", 
    mem="10000",
    species_thresh=parameters.ix['species_number_thresh','entry']
  threads: 1
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    # print(output_on_scratch)
    cp_to_scratch(input, scratch)
    # first pass to get top n species
    shell("touch {output_on_scratch}")
    for i in range(len(input)):
      # print(input_on_scratch[i])
      if file_empty([input_on_scratch[i]]):
        # assert(file_empty([input_on_scratch[i]])),"Input file is empty."
        sample_name = input_on_scratch[i].split('.')[-2]
        filename = input_on_scratch[i]
        shell("sort {filename} > {scratch}/sorted_blast_results.txt")
        with open(scratch+"/sorted_blast_results.txt",'r') as finput:
          last_qseqid = ''
          last_sseqid = ''
          species = []
          for line in finput:
            # must split line with tab because there are spaces in fields
            if line.split('\t')[9] == "Bacteria":
              if line.split('\t')[0] == last_qseqid:
                if line.split('\t')[1] != last_sseqid:
                  print('Warning: '+last_qseqid+' got matched to '+last_sseqid+' and '+line.split()[1]+'\n') 
              else:
                # print(line.split('\t')[7])
                last_qseqid = line.split('\t')[0]
                last_sseqid = line.split('\t')[1]
                species.append(line.split('\t')[7]) # 8th column is the scientific name
        species_cnt = dict([(i,species.count(i)) for i in set(species)])
        # sort species_cnt dictionary by the values
        sorted_keys = sorted(species_cnt, key=species_cnt.get, reverse=True)
        # Append to output file the new subsample
        with open(output_on_scratch[0],'a') as f:
          # The '>' is added in front for later parsing with python script
          t = f.write('>'+sample_name+'\n')
          t = f.write('total\t'+str(len(species))+'\n')
          for k in sorted_keys:
            t = f.write(str(species_cnt[k])+'\t'+k+'\n')
        # shell("echo 'start'")
        # shell("head -n 50 {output_on_scratch[0]}")
        # shell("echo 'end'")
      else:
        sample_name = input_on_scratch[i].split('.')[-2]
        filename = input_on_scratch[i]
        print('Warning: Input file '+filename+' is empty.\n')
        with open(output_on_scratch[0],'a') as f:
          # The '>' is added in front for later parsing with python script
          t = f.write('>'+sample_name+'\n')
          t = f.write('total\t0\n')
    # second path to organize data for plotting
    shell("source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/ &&\
      python {code_dir}/snakehelper_subsample_topspecies.py {output_on_scratch}")
    assert(file_empty([output_on_scratch[0]])),"Output is empty"
    cp_from_scratch(output, scratch)

"""
