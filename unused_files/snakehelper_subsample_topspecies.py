#!/usr/bin/python

"""
change the contig indices to number corresponding to row and col.
This is for matlab plotting so use 1 indexed values

Edit History:
2015.05.25 Created
            
"""

import os, sys, string, glob
import pandas as pd

usage = "usage: snakehelper_subsample_topspecies.py output"

# Check input arguments
if len(sys.argv) < 2:
    print usage
    sys.exit(2)

output = sys.argv[1]

# Read in all lists
lines = []
with open(output, 'r') as f:
  for l in f:
    # print(l)
    lines.append(l)

# Process the lines. Must use tab or new line as delimiter
subsamples = [x[1:].split()[0] for x in lines if '>' in x]
total = [int(x.split()[1]) for x in lines if 'total' in x]
# print(lines)
species = [x.split()[1] for x in lines if '>' not in x and 'total' not in x]
species = list(set(species))
abundance = [[0 for i in range(len(subsamples))] for j in range(len(species))]
for l in lines:
  if '>' in l:
    sample = l[1:].split()[0]
  # need to handle total is 0 case, but it should be fine
  elif 'total' not in l:
    abundance[species.index(l.split()[1])][subsamples.index(sample)] = float(l.split()[0])/float(total[subsamples.index(sample)])

# Write to output
with open(output, 'w') as f:
  f.write('\t'.join(['species']+subsamples)+'\n')
  for i in range(len(species)):
    f.write('\t'.join([species[i]]+[str(x) for x in abundance[i]])+'\n')

print "Organization completed."

