#!/usr/bin/python

"""
change the contig indices to number corresponding to row and col.
This is for matlab plotting so use 1 indexed values
We want each contig to blast against many other contigs so we
do not take only the best match. Instead we take all matches

Edit History:
2015.05.25 Created
2015.05.28 Changed file open synatx to with open() as f
            
"""

import os, sys, string, glob

usage = "usage: snakehelper_contig_similarity.py supercontigname subsamplecontigname blast_report_dist output"

# Check input arguments
if len(sys.argv) < 5:
    print usage
    sys.exit(2)

super_contig_name_file = sys.argv[1]
subsample_name_file = sys.argv[2]
blast_report_file = sys.argv[3]
output = sys.argv[4]

# Read in all lists

super_contig_name = []
with open(super_contig_name_file, 'r') as f:
  for l in f:
    super_contig_name.append(l.split()[0])

subsample_name = []
with open(subsample_name_file, 'r') as f:
  for l in f:
    subsample_name.append(l.split()[0])

a = [] # these are terrible names
b = []
c = []
with open(blast_report_file, 'r') as f:
  for l in f:
    # reading all as strings
    a.append(l.split()[0]) # this is query id
    b.append(l.split()[1]) # this is reference id
    c.append(l.split()[2]) # this is bitscore (the 3rd column)

# Update blast results in a new matrix MATLAB index
assert(len(a) == len(b))
assert(len(b) == len(c))
a_hat = [super_contig_name.index(x)+1 for x in a]
b_hat = [subsample_name.index(x)+1 for x in b]

# handle duplicate situations
blast_dict = {}
for i in range(len(a_hat)):
  # doesn't matter about length of match, just record bitscore
  if (a_hat[i], b_hat[i]) in blast_dict.keys():
    blast_dict[(a_hat[i], b_hat[i])] = max(blast_dict[(a_hat[i], b_hat[i])], c[i])
  else:
    blast_dict[(a_hat[i], b_hat[i])] = c[i]

# Write to output
with open(output, 'w') as f:
  for k in blast_dict.keys():
    f.write('%s\t%s\t%s\n' %(k[0],k[1],blast_dict[k]))

print "File creation completed."

