#!/use/bin/python

# This file moves all the corrected fastq reads from the 
# spade_output_ILxxxxxx/corrected folder into the ILxxxxxx
# root folder. This script moves the files (not copy) so it
# is destructive. It only looks for the .gz versions.
# This script is a patch for transitioning snakemake files
# This script is run from the snakeresults/sample/ folder 

import os, subprocess, sys

usage = "usage: move_corrected_fastq_files.py biosample_result_dir"

if len(sys.argv) < 2:
  print usage
  sys.exit(2)

resultdir = sys.argv[1]

alldir = os.listdir(resultdir)
os.chdir(resultdir)

subsample_name = [x for x in alldir if x[0:2]=='IL']

for folder in subsample_name:

  os.chdir(folder+'/spade_output_'+folder+'/corrected/')

  a=subprocess.call(['mv','ClustPair1.'+folder+'.00.0_0.cor.fastq.gz','../../P1_corrected.'+folder+'.fastq.gz'])
  if a==0:
    print 'file P1_corrected.'+folder+'.fastq.gz moved'
  else:
    print 'file P1_corrected.'+folder+'.fastq.gz not moved'

  a=subprocess.call(['mv','ClustPair2.'+folder+'.00.0_0.cor.fastq.gz','../../P2_corrected.'+folder+'.fastq.gz'])
  if a==0:
    print 'file P2_corrected.'+folder+'.fastq.gz moved'
  else:
    print 'file P2_corrected.'+folder+'.fastq.gz not moved'  

  a=subprocess.call(['mv','ClustPair_unpaired.00.0_0.cor.fastq.gz','../../S_corrected.'+folder+'.fastq.gz'])
  if a==0:
    print 'file S_corrected.'+folder+'.fastq.gz moved'
  else:
    print 'file S_corrected.'+folder+'.fastq.gz not moved'

  os.chdir('../../..')


