#!/bin/bash

################################################################################
# This file copies contigs assembled from the Sherlock cluster.
# It copies the spade output folder from the biosample directories
# on the vm to /local10G/brianyu/snakemake_results/biosample/Combined_Analysis/
# Then it brings the contigs and scaffolds out of the subdirectories.
#################################################################################

# ie. /datastore/brianyu/2015.03.11_SygnisTruePrime/results/spade_output_Sygnis_TruePrime_3species_supercontigs
source_folder=$1

# ie. /local10G/brianyu/snakemake_results/Sygnis_TruePrime_3species/Combined_Analysis/spade_output_Sygnis_TruePrime_3species_supercontigs
destination=$2

parameter_file=$3


biosample=$( head -n 2 $parameter_file | tail -n 1 | cut -f 2 )

rsync -avrP $source_folder $destination

# The corrected reads are already in the Combined_Analysis folder
cd $destination
cp contigs.fasta "../contigs."$biosample".fasta"
cp scaffolds.fasta "../scaffolds."$biosample".fasta"
cd ..


