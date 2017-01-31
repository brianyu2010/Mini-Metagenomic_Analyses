#!/usr/bin/python

"""
This function is used to construct a subsamples.txt file for snakefiles to process
This function should be stored with snakefiles and called with full path to folders
The output should be stored in the code_analysis folder in each biosample
Could later add in lines at the end in order to use subsamples from other runs.

Edit History:
2015.05.20 Created
2015.05.23 Added read count threshold
            
"""

import os, sys, string, glob

usage = "usage: generate_snakemake_seed.py biosample_dir subsamples_file, read_thresh"

# Check input arguments
if len(sys.argv) < 4:
    print usage
    sys.exit(2)

biosample_dir = sys.argv[1]
subsamples_file = sys.argv[2]
read_thresh = int(sys.argv[3])

# Change to folder with all sample folders
print biosample_dir
os.chdir(biosample_dir)

# Read in a list of sequencing run folders
seqruns = os.listdir('.')
seqruns = [seqruns[i] for i in range(len(seqruns)) if seqruns[i][0:6].isdigit()]

# make a python dictionary to keep track of read counts
read_count = {}
seedfile_lines = []

# For each sequencing run
for seqrun in seqruns:

  print seqrun
  
  # Go into the sequencing run folder
  os.chdir(seqrun)
  
  # read in a list of subsamples
  subsamplenames = os.listdir(".")
  # This line works with both Miseq and Nextseq naming conventions
  subsampleIDs = [i.split('__')[-1] for i in subsamplenames]
  print len(subsamplenames)
  print len(subsampleIDs)
  assert len(subsamplenames) == len(subsampleIDs)
  
  # For each subsample
  for i in range(len(subsamplenames)):
    
    # enter a line corresponding to that subsample
    # subsample, biosample, sequencing run
    row = '%s\t%s\t%s\t%s\n' %(subsampleIDs[i],biosample_dir,seqrun,subsamplenames[i])
    seedfile_lines.append(row)

    # check how many lines there are and update
    os.chdir(subsamplenames[i])
    r = os.popen("zcat *R1_001.fastq.gz | wc -l").read().split()[0]
    r = int(r)/4
    print subsampleIDs[i]+" contains "+str(r)+" reads"
    if subsampleIDs[i] in read_count.keys():
      read_count[subsampleIDs[i]] += r
    else:
      read_count[subsampleIDs[i]] = r
    os.chdir('..')
    
  os.chdir('..') # goes back to the biosample level

# threshold read counts, lines
read_count = dict((k, read_count[k]) for k in read_count.keys() if read_count[k] >= read_thresh)
seedfile_lines = [seedfile_lines[i] for i in range(len(seedfile_lines)) if seedfile_lines[i].split()[0] in read_count.keys()]

# Open and prepare output file
outfile = open(subsamples_file,'wb')
header = "%s\t%s\t%s\t%s\n" %('subsampleID','biosample','sequencingrun','subsamplename')
outfile.write(header)
for row in seedfile_lines:
  outfile.write(row)
outfile.close()

outfile = open(biosample_dir + "/results/snakemake_reads.cnt", 'wb')
for k in read_count.keys():
  outfile.write("%s\t%d\n" %(k, read_count[k]))
outfile.close()

print "Snakemake Subsamples.txt Preparation Completed"



